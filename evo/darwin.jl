using Darwin
import Darwin: uniform_mutation

struct NeurodevoInd <: Individual
    genes::Array{Array{Float64}}
    fitness::Array{Float64}
end

function NeurodevoInd(cfg::Dict)
    neuro_cfg = Config("cfg/evo.yaml")
    chromos = make_chromos!(neuro_cfg)
    genes = Array{Array{Float64}}(undef, 10)
    for i in eachindex(chromos)
        genes[i] = chromos[i].genes
    end
    fitness = -Inf*ones(cfg["n_fitness"])
    NeurodevoInd(genes, fitness)
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
