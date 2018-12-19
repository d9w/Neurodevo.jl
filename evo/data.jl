using Neurodevo
using Random

include("cgp.jl")
include("darwin.jl")
include("datasets.jl")

function classify(m::Model, X::Array{Float64}, Y::Array{Int64})
    total_time = 0
    total_memory = 0
    fit = 0.0
    for epoch in 1:m.cfg["n_epochs"]
        labels = zeros(Int64, length(Y))
        for i in eachindex(Y)
            outputs, t, bytes, gctime, memallocs = @timed Neurodevo.step!(m, X[:, i])
            labels[i] = argmax(outputs)
            if labels[i] == Y[i]
                Neurodevo.reward!(m, [1.0])
            end
            total_time += t
            total_memory += bytes
            if (total_time >= m.cfg["time_max"] ||
                total_memory >= m.cfg["memory_max"])
                return [fit - 1.0]
            end
        end
        fit = sum(labels .== Y) / length(Y)
    end
    [fit]
end

function data_evaluation(ind::NeurodevoInd)
    cfg = cgp_cfg(Config("cfg/evo.yaml"), ind.genes[1])
    c = cgp_controller(cfg, ind.genes[2:end])
    m = Model(cfg, c)

    X, Y = ind.func(ind.seed)
    nin = size(X, 1)
    nout = length(unique(Y))
    if cfg["init_method"] == 0
        random_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    else
        layered_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    end
    develop!(m)
    classify(m, X, Y)
end

function NeurodevoInd(genes::Array{Array{Float64}}, fitness::Array{Float64})
    NeurodevoInd(genes, fitness, get_iris, 1)
end

function generation(e::Evolution)
    fits = [get_iris, get_diabetes, get_glass]
    for i in e.population
        i.seed = floor(Int64, (e.gen - 1) / 10)
        func = mod(floor(Int64, (e.gen - 1) / e.cfg["goal_gen"]), length(fits)) + 1
        i.func = fits[func]
    end
end
