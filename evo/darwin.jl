using Darwin
import Darwin: uniform_mutation
import Base.isless
import JSON.json

mutable struct NeurodevoInd <: Individual
    genes::Array{Array{Float64}}
    fitness::Array{Float64}
    func::Function
    seed::Int64
end

function NeurodevoInd(cfg::Dict)
    neuro_cfg = Config("cfg/evo.yaml")
    chromos = make_chromos!(neuro_cfg)
    genes = Array{Array{Float64}}(undef, 0)
    push!(genes, rand(5))
    for c in chromos
        push!(genes, c.genes)
    end
    fitness = -Inf*ones(cfg["n_fitness"])
    NeurodevoInd(genes, fitness)
end

function isless(i1::NeurodevoInd, i2::NeurodevoInd)
    minimum(i1.fitness) < minimum(i2.fitness)
end

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
