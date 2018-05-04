using Gadfly
using Distributions
using Colors
using DataFrames
using Query

include("graph_utils.jl")
include("experiments/seeds.jl")

Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=14pt, point_size=1.2mm,
                        default_color=colorant"#000000"))

colors = [colorant"#1f78b4", colorant"#33a02c", colorant"#e31a1c",
          colorant"#ff7f00", colorant"#6a3d9a", colorant"#a6cee3",
          colorant"#b2df8a", colorant"#fb9a99", colorant"#fdbf6f",
          colorant"#cab2d6"]

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/functions.yaml")

function get_res(log::String)
    res = readtable(log, header=false, separator=' ', names=[
        :epoch, :seed, :Dataset, :Model, :MID, :ARI, :RI, :MI, :HI])
    res[:sid] = string.(res[:Model], "_", res[:MID])
    res
end

function get_bests(res::DataFrame)
    clmeans = by(res, [:Dataset, :Model, :MID], df->DataFrame(arimean=mean(df[:ARI]), aristd=std(df[:ARI])))
end

function plot_acc(res::DataFrame, title::String; filename::String="clustering.pdf",
                  problem::String="iris", ymin=-Inf, ymax=-Inf,
                  ylabel="ARI", key_position=:right)
    cres = @from x in res begin
        @where (x.Dataset == problem)
        @select x; @collect DataFrame
    end
    if ymin == -Inf
        ymin=minimum(cres[:ARI])
    end
    if ymax == -Inf
        ymax=maximum(cres[:ARI])
    end
    cols = [colors[1], colors[2], colors[3]]
    plt = plot(cres, x=:sid, y=:ARI, color=:Model,
               Geom.boxplot,
               Guide.ylabel(ylabel),
               Guide.xlabel("Model"),
               Guide.xticks(label=false),
               Guide.title(title),
               Coord.cartesian(ymin=ymin, ymax=ymax),
               Scale.color_discrete_manual(cols...),
               style(key_position=key_position))

    draw(PDF(filename, 10inch, 4inch), plt)
    plt
end

function pvalues(res::DataFrame, problem::String="iris", model::String="evo0")
    cres = @from x in res begin
        @where (x.Dataset == problem)
        @select x; @collect DataFrame
    end
    lif = @from x in cres begin
        @where (x.Model == "LIF")
        @select x; @collect DataFrame
    end
    pvalues = zeros(10)
    model = @from x in cres begin
        @where (x.Model == model)
        @select x; @collect DataFrame
    end
    for id in 1:10
        x = [lif[:ARI][1:10]; lif[:ARI][1:10]]
        y = model[:ARI][model[:MID] .== id]
        println(id)
        println(size(x))
        println(size(y))
        if size(x) == size(y)
            pvalues[id] = pvalue(OneSampleTTest(Array{Float64}(x), Array{Float64}(y)))
        end
    end
    pvalues
end

function plot_seed10(res::DataFrame, title::String; filename::String="clustering.pdf",
                     ymin=-Inf, ymax=-Inf,
                     ylabel="ARI", key_position=:right)
    plts = Array{Gadfly.Plot}(3)
    problems = unique(res[:Dataset])
    for pid in eachindex(problems)
        problem = problems[pid]
        cres = @from x in res begin
            @where ((x.seed == 10) && (x.Dataset == problem))
            @select x; @collect DataFrame
        end
        cols = [colors[1], colors[2], colors[3]]
        key_position=:none
        if pid == 3
            key_position=:right
        end
        ylabel = nothing
        if pid == 1
            ylabel = "ARI"
        end
        plt = plot(cres, x=:Model, y=:ARI, color=:Model,
                   Geom.boxplot,
                   Guide.ylabel(ylabel),
                   Guide.xlabel(problem),
                   Scale.color_discrete_manual(cols...),
                   style(key_position=:none))
        plts[pid] = plt
    end
    plt = hstack(plts...)

    draw(PDF(filename, 10inch, 4inch), plt)
    plt
end

function plot_graphs(efile::String, edir::String)
    g = lif_graph(-0.65, 0.3, out3=-2, out4=-1)
    c = to_chromo(g)
    chromo_draw(c, "lif.pdf")
    experts = readdlm(efile, ',')
    for exp in 1:size(experts, 1)
        c = PCGPChromo(experts[exp, :], 5, 5)
        chromo_draw(c, string(edir, "/", exp, ".pdf"))
    end
end
