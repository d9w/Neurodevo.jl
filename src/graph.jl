using LightGraphs

function graph_update!(model::Model)
  id = 1
  id_mapping = Dict{Int64, Int64}()
  for (a,s2) in model.synapse
    for s in [model.cells[a].p_id;s2]
      if ~haskey(id_mapping, s)
        id_mapping[s] = id
        id += 1
      end
    end
  end
  model.synapse_graph = DiGraph(id)
  for (a,s) in model.synapse
    add_edge!(model.synapse_graph, id_mapping[model.cells[a].p_id], id_mapping[s])
  end
end

function graph_eval(model::Model)
  cc = global_clustering_coefficient(model.synapse_graph)
  communities = strongly_connected_components(model.synapse_graph)
  divis = ones(Int64, nv(model.synapse_graph))
  for comm=1:length(communities)
    for vert in comm
      divis[vert] = comm
    end
  end
  modul = modularity(Graph(model.synapse_graph), divis)
  dens = LightGraphs.density(model.synapse_graph)
  (nv(model.synapse_graph), ne(model.synapse_graph), cc, length(communities), modul, dens)
end
