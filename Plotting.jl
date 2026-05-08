module Plotting

using CairoMakie
using LaTeXStrings
using ColorSchemes

function color_from_value(x;
    colorrange::Tuple{<:Real,<:Real} = (0.0, 1.0),
    colormap = :viridis,                 # Symbol or ColorScheme
    colorscale = identity,               # e.g., log10
    alpha::Real = 1.0
)
    vmin, vmax = colorrange
    fx, fmin, fmax = colorscale(x), colorscale(vmin), colorscale(vmax)
    t = clamp((fx - fmin) / (fmax - fmin), 0.0, 1.0)

    cs = colormap isa Symbol ? getfield(ColorSchemes, colormap) : colormap  # ColorScheme
    c  = ColorSchemes.get(cs, t)  # RGB
    return RGBAf(c, alpha)
end


export inch, cm, pt, column_width, colors_list, markers_list, get_axis, color_from_value

const inch::Float64 = 96
const cm::Float64 = inch / 2.54
const pt::Float64 = 4 / 3
const column_width::Float64 = 8.6cm
const colors_list = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :gray, :black, :teal, :lime, :gold, :coral, :navy]

const markers_list = [
    :hexagon,
    :circle,
    :pentagon,
    :diamond,
    :star4,
    :vline,
    :cross,
    :xcross,
    :rect,
    :ltriangle,
    :dtriangle,
    :utriangle,
    :star5,
    :star8,
    :star6,
    :rtriangle,
    :octagon,
    :x]

function get_axis(f, row_iter, col_iter, xlabel, ylabel, with_x_labels, with_y_labels, xlimits, ylimits; num_y__minor_ticks::Int64=4, num_x_minor_ticks::Int64=4, title="", titlefontsize::Float64=12pt)
    ax = Axis(f[row_iter, col_iter], xlabel=xlabel, ylabel=ylabel, xlabelsize=12pt, ylabelsize=12pt, xticklabelsize=9pt, yticklabelsize=9pt, spinewidth=2.0pt, xlabelfont=:bold, ylabelfont=:bold, xminorticksvisible=true, yminorticksvisible=true, ylabelvisible=with_y_labels, xlabelvisible=with_x_labels, xticksvisible=true, yticksvisible=true, xticklabelsvisible=with_x_labels, bottomspinevisible=true, xscale=identity, yscale=identity, xtickalign=1.0, xminortickalign=1.0, xtickwidth=2.0pt, xticksize=10.0pt, ytickalign=1.0, yminortickalign=1.0, ytickwidth=2.0pt, yticksize=9.0pt, limits=(xlimits[1], xlimits[2], ylimits[1], ylimits[2]), xminorticks=IntervalsBetween(num_x_minor_ticks), yminorticks=IntervalsBetween(num_y__minor_ticks), xminorticksize=5.0pt, xminortickwidth=2.0pt, yminorticksize=5.0pt, yminortickwidth=2.0pt, title=title, titlesize=titlefontsize)
    return ax
end

end # module Plotting
