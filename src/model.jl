using LightGraphs

include("neuron.jl")

immutable Constants
  morphogen_decay::Float64
  step_init::Int64
  step_total::Int64
  change_stop::Float64
  axon_max::Int64
end

type Model
  dims::Vector{Int64}
  morphogens::Array{Float64}
  neurons::Vector{Neuron}
  neuron_map::Array{Int64}
  synapses::DiGraph
end

function Model(dims::Vector{Int64}=[10,10,10], N::Int64=100)
  morphogens = zeros(Float64, [dims[:];3]...)
  neuron_map = zeros(Float64, dims...)
  neurons = Vector{Neuron}(N)
  for i=1:N
    neurons[i] = Neuron(dims)
    neuron_map[neurons[i].position...] = i
  end
  Model(dims, morphogens, neurons, neuron_map, DiGraph())
end

function step!(model::Model)
  # update morphogens
  for neuron in model.neurons
    ms = model.morphogens[[neuron.position;:]...][:]
    model.morphogens[[neuron.position;:]...] = emission(neuron, ms) + constants.morphogen_decay*ms
  end

  # update neurons
  for neuron in model.neurons
    ms = model.morphogens[[neuron.position;:]...][:]
    action!(neuron, ms, model.dims)
  end

  # update graph
  model.synapses = DiGraph(length(model.neurons))
  for i=1:length(model.neurons)
    neuron = model.neurons[i]
    for axon in neuron.axons
      n = model.neuron_map[axon.position...]
      if n >= 1 && n <= length(model.neurons) && n != i
        add_edge!(model.synapses, i, n)
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
