using Darwin
using RoboGrid
import Darwin: uniform_mutation
import RoboGrid: water_meta, cross_meta
import Base.isless
import JSON.json

include("cgp.jl")
include("darwin.jl")

function NeurodevoInd(genes::Array{Array{Float64}}, fitness::Array{Float64})
    NeurodevoInd(genes, fitness, RoboGrid.coverage_fitness, 1)
end

function generation(e::Evolution)
    fits = [RoboGrid.coverage_fitness, RoboGrid.water_memorize_fitness,
            RoboGrid.cross_memorize_fitness, RoboGrid.cross_strategy_fitness]
    for i in e.population
        i.seed = floor(Int64, (e.gen - 1) / 10)
        func = mod(floor(Int64, (e.gen - 1) / e.cfg["goal_gen"]), length(fits)) + 1
        i.func = fits[func]
    end
end

function robo_controller(m::Model, inputs::Array{Float64})
    Neurodevo.reward!(m, [inputs[1]])
    Neurodevo.step!(m, inputs[2:end])
end

function robo_eval(ind::Individual)
    cfg = cgp_cfg(Config("cfg/evo.yaml"), ind.genes[1])
    c = cgp_controller(cfg, ind.genes[2:end])
    m = Model(cfg, c)
    nin = 9
    nout = length(RoboGrid.ACTIONS)
    if cfg["init_method"] == 0
        random_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    else
        layered_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    end
    develop!(m)
    f(x::Array{Float64}) = robo_controller(m, x)
    [ind.func(f; seed=ind.seed)]
end
