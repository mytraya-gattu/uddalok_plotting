import Pkg;
Pkg.activate("/Users/aragorn/Documents/code/uddalok_plotting/");

if !isdefined(@__MODULE__, :Plotting)
    include("Plotting.jl");
    using .Plotting;
end

using CairoMakie
using LaTeXStrings
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

"""
    discover_U_mu_tuples(results_dir=RESULTS_DIR)

Scan `results_dir` for `.mat` result files and return the available `(U, mu)`
tuples, sorted first by `U` and then by `mu`.
"""
function discover_U_mu_tuples(results_dir::AbstractString=RESULTS_DIR)
    return [(result.U, result.mu) for result in discover_U_mu_results(results_dir)]
end

"""
    cell_edges(values)

Infer heatmap cell edges from sorted grid-center values.
"""
function cell_edges(values::AbstractVector{<:Real})
    length(values) >= 2 || error("Need at least two grid values to infer cell edges.")

    edges = Vector{Float64}(undef, length(values) + 1)
    for i in 2:length(values)
        edges[i] = (values[i - 1] + values[i]) / 2
    end

    edges[1] = values[1] - (values[2] - values[1]) / 2
    edges[end] = values[end] + (values[end] - values[end - 1]) / 2
    return edges
end

"""
    shade_odd_chern_boxes!(ax, U_values, mu_values, chern_grid, colorrange, colormap)

Shade every finite odd-Chern cell in the `U, mu` heatmap using hatch fill.
"""
function shade_odd_chern_boxes!(ax, U_values, mu_values, chern_grid, colorrange, colormap)
    U_edges = cell_edges(U_values)
    mu_edges = cell_edges(mu_values)
    hatch = Vec2f(1, 1)

    for i in eachindex(U_values), j in eachindex(mu_values)
        C = chern_grid[i, j]
        isfinite(C) && isodd(round(Int, C)) || continue

        xlo, xhi = U_edges[i], U_edges[i + 1]
        ylo, yhi = mu_edges[j], mu_edges[j + 1]
        corners = Point2f[
            Point2f(xlo, ylo),
            Point2f(xhi, ylo),
            Point2f(xhi, yhi),
            Point2f(xlo, yhi),
        ]

        hatch_color = color_from_value(C; colorrange = colorrange, colormap = colormap, alpha = 0.95)
        pattern = Makie.LinePattern(;
            direction = hatch,
            width = 5,
            tilesize = (20, 20),
            linecolor = hatch_color,
            backgroundcolor = (:white, 0.2))

        poly!(ax, corners;
            color = pattern,
            strokewidth = 1.0pt, strokecolor=:red)
    end

    return ax
end

let
    U_mu_results = discover_U_mu_results()
    @info "Discovered $(length(U_mu_results)) result points for Chern-number map" results_dir=RESULTS_DIR

    U_values = sort!(unique(result.U for result in U_mu_results))
    mu_values = sort!(unique(result.mu for result in U_mu_results))

    U_index = Dict(U => i for (i, U) in enumerate(U_values))
    mu_index = Dict(mu => i for (i, mu) in enumerate(mu_values))
    chern_grid = fill(NaN, length(U_values), length(mu_values))

    for result in U_mu_results
        data = matread(result.path)
        data["delta_max"] < 1e-4 && continue

        chern_grid[U_index[result.U], mu_index[result.mu]] = round(Int, data["C_BdG"])
    end

    chern_numbers = sort!(unique(Int.(filter(isfinite, vec(chern_grid)))))
    isempty(chern_numbers) && error("No finite Chern numbers found after delta_max filtering.")
    chern_color_limits = (minimum(chern_numbers) - 0.5, maximum(chern_numbers) + 0.5)

    fig = Figure(size = (column_width, column_width), figure_padding = (2.0pt, 2.0pt, 2.0pt, 2.0pt))
    ax = get_axis(fig, 1, 1, L"U / t", L"\mu / t", true, true, (nothing, nothing), (nothing, nothing))

    ax.xlabelsize = 12.0pt
    ax.ylabelsize = 12.0pt

    colormap = :viridis
    heatmap!(ax, U_values, mu_values, chern_grid;
        colormap = colormap,
        colorrange = chern_color_limits,
        nan_color = (:grey, 0.35))
    shade_odd_chern_boxes!(ax, U_values, mu_values, chern_grid, chern_color_limits, colormap)

    Colorbar(fig[begin:end, 2],
        colormap = colormap,
        limits = chern_color_limits,
        width = 8pt)
    Label(fig[begin:end, 2, Right()],
        L"C_{\mathrm{BdG}}",
        fontsize = 12.0pt,
        color = :black,
        padding = (2.0pt, 0, 0pt, 0))
    colgap!(fig.layout, 1, 2.0pt)

    hidedecorations!(ax, label = false, ticks = false, ticklabels = false)
    fig
end
