import Pkg;
Pkg.activate("/Users/aragorn/Documents/code/uddalok_plotting/");

include("Plotting.jl");
using .Plotting;

using CairoMakie
using LaTeXStrings
using JLD2
using MAT

const RESULTS_DIR = joinpath(@__DIR__, "results_Wdis_0.05_dis_1_8_8_dense_c")
const U_MU_PATTERN = r"(?:^|_)U=([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)_mu=([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"

"""
    parse_U_mu_tuple(path)

Parse a result filename containing `U=<value>_mu=<value>` and return
`(U, mu)`. Returns `nothing` if the filename does not match the expected
result-file convention.
"""
function parse_U_mu_tuple(path::AbstractString)
    matched = match(U_MU_PATTERN, basename(path))
    isnothing(matched) && return nothing
    return (parse(Float64, matched.captures[1]), parse(Float64, matched.captures[2]))
end

"""
    discover_U_mu_results(results_dir=RESULTS_DIR)

Scan `results_dir` for `.mat` result files and return named tuples containing
the parsed `U`, `mu`, and exact file `path`, sorted first by `U` and then by
`mu`.
"""
function discover_U_mu_results(results_dir::AbstractString=RESULTS_DIR)
    isdir(results_dir) || error("Results directory does not exist: $results_dir")

    results = Dict{Tuple{Float64, Float64}, String}()
    for path in readdir(results_dir; join=true)
        isfile(path) || continue
        endswith(path, ".mat") || continue

        parsed = parse_U_mu_tuple(path)
        isnothing(parsed) && continue
        results[parsed] = path
    end

    pairs = sort!(collect(keys(results)), by = pair -> (pair[1], pair[2]))
    return [(U = pair[1], mu = pair[2], path = results[pair]) for pair in pairs]
end

let

    U_mu_results = discover_U_mu_results()
    @info "Discovered $(length(U_mu_results)) result points for parity contours" results_dir=RESULTS_DIR

    U_values = sort!(unique(result.U for result in U_mu_results))
    μ_values = sort!(unique(result.mu for result in U_mu_results))

    U_index = Dict(U => i for (i, U) in enumerate(U_values))
    μ_index = Dict(μ => i for (i, μ) in enumerate(μ_values))
    parity_grid = fill(NaN, length(U_values), length(μ_values))

    for result in U_mu_results
        data = matread(result.path)
        Δ = data["delta_max"]
        Δ < 1e-4 && continue

        C = round(Int, data["C_BdG"])
        parity_grid[U_index[result.U], μ_index[result.mu]] = isodd(C) ? 1.0 : 0.0
    end

    fig = Figure(size=(column_width, column_width), figure_padding=(2.0pt, 2.0pt, 2.0pt, 2.0pt))
    ax = get_axis(fig, 1, 1, L"U / t", L"\mu / t", true, true, (nothing, nothing), (nothing, nothing))

    ax.xlabelsize = 12.0pt
    ax.ylabelsize = 12.0pt

    contourf!(ax, U_values, μ_values, parity_grid;
        levels = [-0.5, 0.5, 1.5],
        colormap = [RGBAf(to_color(colors_list[2]), 0.35), RGBAf(to_color(colors_list[1]), 0.35)])
    contour!(ax, U_values, μ_values, parity_grid;
        levels = [0.5],
        color = :black,
        linewidth = 1.5pt)

    nan_points = [(result.U, result.mu) for result in U_mu_results if isnan(parity_grid[U_index[result.U], μ_index[result.mu]])]
    scatter!(ax, first.(nan_points), last.(nan_points);
        color = (:grey, 0.5),
        marker = :xcross,
        markersize = 7.0pt,
        strokewidth = 0.75pt)

    odd_marker = MarkerElement(color = RGBAf(to_color(colors_list[1]), 0.35), marker = :rect, markersize = 10.0pt)
    even_marker = MarkerElement(color = RGBAf(to_color(colors_list[2]), 0.35), marker = :rect, markersize = 10.0pt)
    axislegend(ax, [odd_marker, even_marker], ["Odd Chern", "Even Chern"];
        position = (0.05, 0.05),
        padding = (0.0pt, 2.0pt, 0.0pt, 0.0pt),
        rowgap = -2.0pt,
        labelsize = 12pt,
        patchlabelgap = 1.0pt)

    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    fig

end

"""
    discover_U_mu_tuples(results_dir=RESULTS_DIR)

Scan `results_dir` for `.mat` result files and return the available `(U, mu)`
tuples, sorted first by `U` and then by `mu`.
"""
function discover_U_mu_tuples(results_dir::AbstractString=RESULTS_DIR)
    return [(result.U, result.mu) for result in discover_U_mu_results(results_dir)]
end

let 

    U_mu_results = discover_U_mu_results()
    U_mu_tuples = [(result.U, result.mu) for result in U_mu_results]
    @info "Discovered $(length(U_mu_tuples)) result points" results_dir=RESULTS_DIR

    fig = Figure(size=(column_width, column_width), figure_padding=(2.0pt, 2.0pt, 2.0pt, 2.0pt))

    ax = get_axis(fig, 1, 1, L"U / t", L"\mu / t", true, true, (nothing, nothing), (nothing, nothing))

    ax.xlabelsize = 12.0pt
    ax.ylabelsize = 12.0pt


    chern_values = Vector{Float64}(undef, length(U_mu_results))
    

    results = []
    for (i, result) in enumerate(U_mu_results)
        data = matread(result.path)
        chern_values[i] = data["C_BdG"]
        Δ = data["delta_max"]
        if Δ < 1e-4
            chern_values[i] = NaN
        else
            chern_values[i] = round(Int, data["C_BdG"])
        end
        U = result.U
        μ = result.mu
        push!(results, (U, μ, chern_values[i]))
    end

    cmap_val = :viridis

    odd_values = filter(x->isfinite(x[3]) && isodd(x[3]), results)
    even_values = filter(x->isfinite(x[3]) && !isodd(x[3]), results)
    nan_values = filter(x->!isfinite(x[3]), results)

    filter!(x->!(isnan(x)), chern_values)
    helper(c) = isnan(c) ? :grey : color_from_value(c, colorrange=(minimum(chern_values), maximum(chern_values)), colormap=cmap_val)

    scatter!(ax, first.(odd_values),map(x->x[2], odd_values), color=helper.(map(x->x[3], odd_values)), marker=markers_list[1], markersize=10.0pt, strokecolor=:black, strokewidth=1.0pt, label="Odd Chern", colormap=cmap_val)
    scatter!(ax, first.(even_values), map(x->x[2], even_values), color=helper.(map(x->x[3], even_values)), marker=markers_list[2], markersize=10.0pt, strokecolor=:black, strokewidth=1.0pt, label="Even Chern", colormap=cmap_val)

    Colorbar(fig[begin:end, 2], colormap=cmap_val, limits=(minimum(chern_values), maximum(chern_values)), width=8pt)

    Label(fig[begin:end, 2, Right()], L"C_{\mathrm{BdG}}", fontsize=12.0pt, color=:black, padding=(2.0pt, 0, 00pt, 0))
    colgap!(fig.layout, 1, 2.0pt)
    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    fig 
    
    
end
