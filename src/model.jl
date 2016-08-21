using Base.Cartesian
using LightGraphs
using Distances
using Interpolations

include("neuron.jl")

immutable Constants
  morphogen_decay::Float64
  step_init::Int64
  step_total::Int64
  change_stop::Float64
  axon_max::Int64
  neuron_size::Float64
  min_dist::Float64
end

type Model
  dims::Vector{Int64}
  morphogens::Array{Float64}
  grid::Array{Int64}
  neurons::Vector{Neuron}
  synapses::DiGraph
end

function Model(dims::Vector{Int64}=[10,10,10], N::Int64=100)
  morphogens = zeros(Float64, [dims[:];3]...)
  grid = zeros(Int64, *(dims...), length(dims))
  i = 1
  @nloops length(dims) j morphogens begin
    grid[i] = collect(@ntuple length(dims) j)
    i += 1
  end
  neurons = Vector{Neuron}(N)
  for i=1:N
    neurons[i] = Neuron(dims)
  end
  Model(dims, morphogens, neurons, DiGraph())
end

function step!(model::Model)
  # update morphogens
  model.morphogens *= constants.morphogen_decay
  for i in eachindex(model.grid[:,1])
    pos = grid[i,:]
    for neuron in model.neurons
      morphogens[pos...;1] += 1/max(constants.min_dist,evaluate(Euclidean(), pos, neuron.position))
      for axon in neuron.axons
        morphogens[pos...;2] += 1/max(constants.min_dist,evaluate(Euclidean(), pos, axon.position))
      end
    end
    # this can be pre-computed
    morphogens[pos...;3] += 1/max(constants.min_dist,evaluate(Euclidean(), pos, dims/2))
  end

  itp = interpolate(morphogens, BSpline(Linear()), OnGrid())

  # update neurons
  for neuron in model.neurons
    morphogens = [itp[[neuron.position[:];i]...] for i=1:3]
    gradients = [gradient(itp, neuron.position... ,i)[1:length(model.dims)] for i=1:3]
    action!(neuron, ms, gradients, model.dims)
  end

  # update graph
  model.synapses = DiGraph(length(model.neurons))
  for i=1:length(model.neurons)
    neuron = model.neurons[i]
    for axon in neuron.axons
      for j=1:length(model.neurons)
        if evaluate(Euclidean(), axon.position, model.neurons[n].position) < constants.neuron_size
          add_edge!(model.synapses, i, n)
        end
      end
    end
  end
end

function evaluate(model::Model)
  cc = global_clustering_coefficient(model.synapses)
  communities = strongly_connected_components(model.synapses)
  divis = ones(Int64, nv(model.synapses))
  for comm=1:length(communities)
    for vert in comm
      divis[vert] = comm
    end
  end
  modul = modularity(Graph(model.synapses), divis)
  dens = density(model.synapses)
  n_axons = 0
  for n in model.neurons
    n_axons += length(n.axons)
  end
  n_axons /= length(model.neurons)
  (nv(model.synapses), ne(model.synapses), n_axons, cc, length(communities), modul, dens)
end
