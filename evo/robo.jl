using Darwin
using RoboGrid
import Darwin: uniform_mutation
import RoboGrid: water_meta, cross_meta
import Base.isless
import JSON.json

include("cgp.jl")

mutable struct NeurodevoInd <: Individual
    genes::Array{Array{Float64}}
    fitness::Array{Float64}
    func::Function
    seed::Int64
end

function NeurodevoInd(genes::Array{Array{Float64}}, fitness::Array{Float64})
    NeurodevoInd(genes, fitness, RoboGrid.coverage_fitness, 1)
end

function NeurodevoInd(cfg::Dict)
    neuro_cfg = Config("cfg/evo.yaml")
    chromos = make_chromos!(neuro_cfg)
    genes = Array{Array{Float64}}(undef, 0)
    push!(genes, rand(5))
    for c in chromos
        push!(genes, c.genes)
    end
    NeurodevoInd(genes, [-Inf])
end

isless(i1::NeurodevoInd, i2::NeurodevoInd) = i1.fitness[1] < i2.fitness[1]

function json(i::NeurodevoInd)
    json(Dict("genes"=>i.genes, "fitness"=>i.fitness))
end

function uniform_mutation(parent::NeurodevoInd; m_rate=0.1)
    genes = Array{Array{Float64}}(undef, 0)
    for gene in parent.genes
        cgene = rand(length(gene))
        inds = rand(length(gene)) .>= m_rate
        cgene[inds] = gene[inds]
        push!(genes, cgene)
    end
    NeurodevoInd(genes, deepcopy(parent.fitness))
end

function generation(e::Evolution)
    fits = [RoboGrid.coverage_fitness, RoboGrid.water_search_fitness,
            RoboGrid.water_memorize_fitness, RoboGrid.cross_search_fitness,
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
