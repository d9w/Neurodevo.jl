using Gadfly
using Colors
using DataFrames
using Query
using CSV
using Statistics
import Cairo, Fontconfig
import ColorSchemes


Gadfly.push_theme(Theme(major_label_font="Helvetica", major_label_font_size=20pt,
                        minor_label_font="Helvetica", minor_label_font_size=12pt,
                        key_title_font="Helvetica", key_title_font_size=16pt,
                        key_label_font="Helvetica", key_label_font_size=14pt,
                        line_width=0.8mm, point_size=1.0mm,
                        highlight_width=0.1mm,
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        default_color=colorant"#000000"))

colors = [colorant"#a6cee3", colorant"#1f78b4", colorant"#b2df8a",
          colorant"#33a02c", colorant"#fb9a99", colorant"#e31a1c"]

function op_color(p::Int64)
    # eval(Meta.parse(string("ColorSchemes.Set1_", max(p, 3))))
    # reshape(ColorSchemes.Paired_8, (2, 4))'[:]
    ColorSchemes.Set1_4
end

function get_evolution(logfile::String; id::String="", goalgen::Int64=10000)
    res = CSV.File(logfile; delim=" ",
                   header=[:timestamp, :darwin, :info, :seed, :generation,
                           :maxfit, :meanfit, :stdfit]) |> DataFrame
    if length(id) > 0
        res[:id] = id
    else
        res[:id] = string(res[:seed][1])
    end
    fits = ["coverage_fitness", "water_search_fitness",
            "water_memorize_fitness", "cross_search_fitness",
            "cross_memorize_fitness", "cross_strategy_fitness"]
    res[:fit] = mod.(floor.(Int64, (res[:generation] .- 1) / goalgen), 4) .+ 1
    res[:fit_id] = string.(res[:id], ", ", res[:fit])
    res
end

function plot_evolution(res::DataFrame; title::String="Evolution",
                        filename="evolution.pdf")
    plt = plot(res, x=:generation, y=:maxfit, color=:id,
               # shape=:id,
               ygroup=:fit,
               Geom.subplot_grid(free_y_axis=true,
               # Geom.point,
               Geom.smooth(smoothing=0.2)),
               Guide.xlabel("Generation"),
               Guide.ylabel("Fitness"),
               # Scale.color_discrete_hue(op_color),
               Guide.title(title))
    draw(PDF(filename, 8inch, 8inch), plt)
end

function plot_all()
    res = get_evolution("results/88778/evolution.log", id="CGP",
                        goalgen=100)
    res = vcat(res, get_evolution("results/88779/evolution.log", id="Nd",
                                  goalgen=100))
    plot_evolution(res; filename="evolution_100.pdf")

    res = get_evolution("results/88794/evolution.log",
                        id="CGP", goalgen=10)
    res = vcat(res, get_evolution("results/88796/evolution.log",
                                  id="Nd", goalgen=10))
    plot_evolution(res; filename="evolution_10.pdf")
end


