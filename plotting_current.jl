import Pkg;
Pkg.activate("/Users/aragorn/Documents/code/uddalok_plotting/");

if !isdefined(@__MODULE__, :Plotting)
    include("Plotting.jl");
    using .Plotting;
end

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
    @info "Discovered $(length(U_mu_results)) result points for parity surface" results_dir=RESULTS_DIR

    U_values = sort!(unique(result.U for result in U_mu_results))
    μ_values = sort!(unique(result.mu for result in U_mu_results))

    U_index = Dict(U => i for (i, U) in enumerate(U_values))
    μ_index = Dict(μ => i for (i, μ) in enumerate(μ_values))
    parity_height = fill(NaN, length(U_values), length(μ_values))
    chern_grid = fill(NaN, length(U_values), length(μ_values))

    for result in U_mu_results
        data = matread(result.path)
        Δ = data["delta_max"]
        Δ < 1e-4 && continue

        C = round(Int, data["C_BdG"])
        parity_height[U_index[result.U], μ_index[result.mu]] = isodd(C) ? 1.0 : -1.0
        chern_grid[U_index[result.U], μ_index[result.mu]] = C
    end

    chern_numbers, chern_colors, chern_to_color_index = interleaved_chern_colormap(chern_grid)
    color_index_grid = chern_color_index_grid(chern_grid, chern_to_color_index)
    chern_color_limits = (0.5, length(chern_numbers) + 0.5)

    fig = Figure(size=(column_width, column_width), figure_padding=(2.0pt, 2.0pt, 2.0pt, 2.0pt))
    ax = Axis3(fig[1, 1],
        xlabel = L"U / t",
        ylabel = L"\mu / t",
        zlabel = L"\mathrm{parity}",
        xlabelsize = 12.0pt,
        ylabelsize = 12.0pt,
        zlabelsize = 12.0pt,
        azimuth = 0.85pi,
        elevation = 0.22pi,
        protrusions = (0.0pt, 0.0pt, 0.0pt, 0.0pt))

    surface!(ax, U_values, μ_values, parity_height;
        color = color_index_grid,
        colormap = chern_colors,
        colorrange = chern_color_limits,
        nan_color = (:grey, 0.0),
        shading = NoShading)
    wireframe!(ax, U_values, μ_values, parity_height;
        color = (:black, 0.25),
        linewidth = 0.5pt)

    ax.zticks = ([-1.0, 1.0], ["even", "odd"])

    Colorbar(fig[1, 2],
        colormap = chern_colors,
        limits = chern_color_limits,
        ticks = (1:length(chern_numbers), string.(chern_numbers)),
        width = 8pt)
    Label(fig[1, 2, Right()], L"C_{\mathrm{BdG}}", fontsize = 12.0pt, color = :black, padding = (2.0pt, 0, 0pt, 0))
    colgap!(fig.layout, 1, 2.0pt)

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

    odd_values = filter(x->isfinite(x[3]) && isodd(x[3]), results)
    even_values = filter(x->isfinite(x[3]) && !isodd(x[3]), results)
    nan_values = filter(x->!isfinite(x[3]), results)

    filter!(x->!(isnan(x)), chern_values)
    chern_numbers, chern_colors, chern_to_color_index = interleaved_chern_colormap(chern_values)
    chern_color_limits = (0.5, length(chern_numbers) + 0.5)
    helper(c) = isnan(c) ? :grey : chern_colors[chern_to_color_index[Int(c)]]

    scatter!(ax, first.(odd_values),map(x->x[2], odd_values), color=helper.(map(x->x[3], odd_values)), marker=markers_list[1], markersize=10.0pt, strokecolor=:black, strokewidth=1.0pt, label="Odd Chern")
    scatter!(ax, first.(even_values), map(x->x[2], even_values), color=helper.(map(x->x[3], even_values)), marker=markers_list[2], markersize=10.0pt, strokecolor=:black, strokewidth=1.0pt, label="Even Chern")

    Colorbar(fig[begin:end, 2], colormap=chern_colors, limits=chern_color_limits, ticks=(1:length(chern_numbers), string.(chern_numbers)), width=8pt)

    Label(fig[begin:end, 2, Right()], L"C_{\mathrm{BdG}}", fontsize=12.0pt, color=:black, padding=(2.0pt, 0, 00pt, 0))
    colgap!(fig.layout, 1, 2.0pt)
    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    fig 
    
    
end

let

    U_mu_results = discover_U_mu_results()
    @info "Discovered $(length(U_mu_results)) result points for parity contours" results_dir=RESULTS_DIR

    U_values = sort!(unique(result.U for result in U_mu_results))
    μ_values = sort!(unique(result.mu for result in U_mu_results))

    U_index = Dict(U => i for (i, U) in enumerate(U_values))
    μ_index = Dict(μ => i for (i, μ) in enumerate(μ_values))
    parity_grid = fill(NaN, length(U_values), length(μ_values))

    chern_values = fill(NaN, length(U_values), length(μ_values))
    for result in U_mu_results
        data = matread(result.path)
        Δ = data["delta_max"]
        C = round(Int, data["C_BdG"])
        if Δ < 1e-4
            parity_grid[U_index[result.U], μ_index[result.mu]] = -1.0
            chern_values[U_index[result.U], μ_index[result.mu]] = NaN
        elseif isodd(C)
            parity_grid[U_index[result.U], μ_index[result.mu]] = 0.0
            chern_values[U_index[result.U], μ_index[result.mu]] = C

        else
            parity_grid[U_index[result.U], μ_index[result.mu]] = 1.0
            chern_values[U_index[result.U], μ_index[result.mu]] = C
        end
        # parity_grid[U_index[result.U], μ_index[result.mu]] = 
        # chern_values[U_index[result.U], μ_index[result.mu]] = C
    end

    fig = Figure(size=(column_width, column_width), figure_padding=(2.0pt, 2.0pt, 2.0pt, 2.0pt))
    ax = get_axis(fig, 1, 1, L"U / t", L"\mu / t", true, true, (nothing, nothing), (nothing, nothing))

    ax.xlabelsize = 12.0pt
    ax.ylabelsize = 12.0pt

    chern_numbers, chern_colors, chern_to_color_index = interleaved_chern_colormap(chern_values)
    color_index_grid = chern_color_index_grid(chern_values, chern_to_color_index)
    chern_color_limits = (0.5, length(chern_numbers) + 0.5)
    # contourf!(ax, U_values, μ_values, parity_grid;
    #     levels = [-0.5, 0.5, 1.5],
    #     colormap = [RGBAf(to_color(colors_list[2]), 0.35), RGBAf(to_color(colors_list[1]), 0.35)])
    heatmap!(ax, U_values, μ_values, color_index_grid;
        colormap = chern_colors,
        colorrange = chern_color_limits,
        nan_color = (:grey, 0.35))
    contour!(ax, U_values, μ_values, parity_grid; levels=3, linewidth=2.0pt, linestyle=:solid, colormap=:RdBu)

    # nan_points = [(result.U, result.mu) for result in U_mu_results if isnan(parity_grid[U_index[result.U], μ_index[result.mu]])]
    # scatter!(ax, first.(nan_points), last.(nan_points);
    #     color = (:grey, 0.5),
    #     marker = :xcross,
    #     markersize = 7.0pt,
    #     strokewidth = 0.75pt)

    # odd_marker = MarkerElement(color = RGBAf(to_color(colors_list[1]), 0.35), marker = :rect, markersize = 10.0pt)
    # even_marker = MarkerElement(color = RGBAf(to_color(colors_list[2]), 0.35), marker = :rect, markersize = 10.0pt)
    # axislegend(ax, [odd_marker, even_marker], ["Odd Chern", "Even Chern"];
    #     position = (0.05, 0.05),
    #     padding = (0.0pt, 2.0pt, 0.0pt, 0.0pt),
    #     rowgap = -2.0pt,
    #     labelsize = 12pt,
    #     patchlabelgap = 1.0pt)

    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    fig

end

let

    U_mu_results = discover_U_mu_results()
    @info "Discovered $(length(U_mu_results)) result points for transition boxes" results_dir=RESULTS_DIR

    U_values = sort!(unique(result.U for result in U_mu_results))
    μ_values = sort!(unique(result.mu for result in U_mu_results))

    U_index = Dict(U => i for (i, U) in enumerate(U_values))
    μ_index = Dict(μ => i for (i, μ) in enumerate(μ_values))

    parity_grid = fill(NaN, length(U_values), length(μ_values))
    chern_values = fill(NaN, length(U_values), length(μ_values))

    for result in U_mu_results
        data = matread(result.path)
        Δ = data["delta_max"]
        Δ < 1e-4 && continue

        C = round(Int, data["C_BdG"])
        parity_grid[U_index[result.U], μ_index[result.mu]] = isodd(C) ? 1.0 : 0.0
        chern_values[U_index[result.U], μ_index[result.mu]] = C
    end

    fig = Figure(size=(column_width, column_width), figure_padding=(2.0pt, 2.0pt, 2.0pt, 2.0pt))
    ax = get_axis(fig, 1, 1, L"U / t", L"\mu / t", true, true, (nothing, nothing), (nothing, nothing))

    ax.xlabelsize = 12.0pt
    ax.ylabelsize = 12.0pt

    chern_numbers, chern_colors, chern_to_color_index = interleaved_chern_colormap(chern_values)
    color_index_grid = chern_color_index_grid(chern_values, chern_to_color_index)
    chern_color_limits = (0.5, length(chern_numbers) + 0.5)
    heatmap!(ax, U_values, μ_values, color_index_grid;
        colormap = chern_colors,
        colorrange = chern_color_limits,
        nan_color = (:grey, 0.35))

    transition_boxes = Tuple{Float64, Float64, Float64, Float64}[]
    for i in 1:(length(U_values) - 1), j in 1:(length(μ_values) - 1)
        cell_parities = [
            parity_grid[i, j],
            parity_grid[i + 1, j],
            parity_grid[i, j + 1],
            parity_grid[i + 1, j + 1],
        ]
        finite_parities = unique(filter(isfinite, cell_parities))
        length(finite_parities) > 1 || continue
        push!(transition_boxes, (U_values[i], U_values[i + 1], μ_values[j], μ_values[j + 1]))
    end
    Colorbar(fig[begin:end, 2],
        colormap = chern_colors,
        limits = chern_color_limits,
        ticks = (1:length(chern_numbers), string.(chern_numbers)),
        width = 8pt)
    Label(fig[begin:end, 2, Right()], L"C_{\mathrm{BdG}}", fontsize=12.0pt, color=:black, padding=(2.0pt, 0, 0pt, 0))
    colgap!(fig.layout, 1, 2.0pt)

    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    fig

end
