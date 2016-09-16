using LightGraphs

function graph_update!(model::Model)
  model.synapses = DiGraph(length(model.neurons))
  for i=1:length(model.neurons)
    neuron = model.neurons[i]
    for axon in neuron.axons
      for j=1:length(model.neurons)
        if evaluate(Euclidean(), axon.position, model.neurons[j].position) < constants.neuron_size
          add_edge!(model.synapses, i, j)
        end
      end
    end
  end
end

function graph_eval(model::Model)
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
