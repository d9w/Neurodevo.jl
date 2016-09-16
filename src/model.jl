using LightGraphs
using Distances
using Grid

include("controller.jl")

type Cell
  p_id::Int64
  id::Int64
  position::Vector{Float64}
  params::Vector{Int64}
  ctype::Int64
end

function Cell(id::Int64)
  position = DIMS[:,1]+rand(size(DIMS,1)).*(DIMS[:,2]-DIMS[:,1])
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(0, id, position, params, ctype)
end

type Model
  maxid::Int64
  morphogens::Array{Float64}
  grid::Array{Int64}
  cells::Dict{Int64, Cell}
  soma_axons::Dict{Int64, Vector{Int64}}
  synapse::Dict{Int64, Int64}
  synapse_weights::Dict{Int64, Float64}
  synapse_graph::DiGraph
end

function Model()
  i = 1
  grid = Array{Float64}(*(DIMS...), N_D)
  for c = Counter(DIMS)
    grid[i,:] = convert(Array{Float64}, c)
    i += 1
  end
  morphogens = zeros(Float64, [DIMS[:];N_MORPHS]...)
  cells = Dict{ASCIIString, Cell}
  for i=1:N_NSC
    cell = Cell(i)
    cells[i] = cell
  end
  Model(N_NSC, morphogens, grid, cells, Dict{Int64, Vector{Int64}}(), Dict{Int64, Int64}(), DiGraph())
end

function apoptosis!(model::Model, cell::Cell)
  if cell.ctype == 3
    delete!(model.soma_axons, cell.id)
    for s in keys(model.synapses)
      if model.synapses[s] == cell.id
        delete!(model.synapses, s)
        delete!(model.synapse_weights, s)
      end
    end
  elseif cell.ctype == 4
    delete!(model.synapses, cell.id)
    delete!(model.synapse_weights, cell.id)
    for s in keys(model.soma_axons)
      for a in model.soma_axons[s]
        if a == cell.id
          pop!(model.soma_axons[s],a)
          break
        end
      end
    end
  end
  delete!(model.cells, cell.id)
end

function update_morphogens!(model::Model)
  # model.morphogens *= MORPHOGEN_DECAY
  for i in eachindex(model.grid[:,1])
    pos = vec(model.grid[i,:])
    for (ckey, cell) in cells
      for m in 1:N_MORPHS
        morph = model.morphogens[[pos;m]...]
        morph += morphogen_diff(m, cell.type, cell.params, cell.position, pos)
        morph = max(0.0, morph)
        model.morphogens[[pos;m]...] = morph
      end
    end
  end
  # morphogen normalization
  # sm = [size(morphs)[1:N_D]...]
  # maxmorphs = vec(maximum(model.morphogens, collect(1:N_D)))
  # for m in eachindex(maxmorphs)
  #   model.morphogens[[[1:sm[i] for i in eachindex(sm)];m]...] /= maxmorphs[m]
  # end
end

function cell_division!(model::Model, itp::Grid.InterpGrid)
  for (ckey,cell) in model.cells
    morphogens = [itp[[cell.position[:];m]...] for m in 1:N_MORPHS]
    if division(morphogens, cell.ctype, cell.parameters)
      model.maxid += 1
      branch_dir = child_branch(morphogens, cell.ctype, cell.parameters)
      ctype = child_type(morphogens, cell.ctype, cell.parameters)
      id = model.maxid
      if ctype == 0
        apoptosis!(model, cell)
      else
        params = child_parameters(morphogens, cell.ctype, cell.parameters, ctype)
        pos = cell.position + child_position(morphogens, cell.ctype, cell.parameters, DIMS)
        pos = min(max(pos, [0.,0.,0.]), DIMS)
        model.cells[id] = Cell(cell.id, id, pos, params, ctype)
        if ctype == 3
          model.soma_axons[id] = []
        elseif cell.ctype == 3 && ctype == 4
          push!(model.soma_axons[cell.id], id)
        end
      end
    end
  end
end

function synapse_update!(model::Model, itp::Grid.InterpGrid)
  for (akey, acell) in model.cells
    if acell.ctype == 4
      if ~haskey(akey, model.synapses)
        for (skey, scell) in model.cells
          if scell.ctype == 3 && bkey != acell.p_id
            if synapse_formation(acell.position, bcell.position, DIMS)
              model.synapses[akey] = skey
              model.synapse_weights[akey] = 0.0
              break
            end
          end
        end
      end
      if haskey(akey, model.synapses)
        soma = model.cells[model.synapse[akey]]
        amorphs = [itp[[acell.position[:];m]...] for m in 1:N_MORPHS]
        smorphs = [itp[[soma.position[:];m]...] for m in 1:N_MORPHS]
        model.synapse_weights[akey] += synapse_weight(smorphs, amorphs, soma.params, acell.params)
        if ~synapse_survival(model.synapse_weights[akey])
          apoptosis!(model, acell)
        end
      end
    end
  end
end

function cell_movement!(model::Model, itp::Grid.InterpGrid)
  for (ckey,cell) in model.cells
    morphs = Vector{Float64}(N_MORPHS)
    grad = Array{Float64}(N_MORPHS, N_D)
    for m = 1:N_MORPHS
      v,g = valgrad(msi, [pos[:];m]...)
      morps[m] = v
      grad[m,:] = g[1:N_D]
    end
    if ~in(ckey, [collect(keys(model.synapses));collect(values(model.synapses))])
      cell.position += cell_movement(morphs, grad, cell.ctype, cell.params, DIMS)
      cell.position = min(max(cell.position, [0.,0.,0.]), DIMS)
    end
  end
end

function step!(model::Model)
  # update morphogens
  update_morphogens!(model)
  itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)

  # division and apoptosis
  cell_division!(model, itp)

  # synapse formation and checking
  synapse_update!(model, itp)

  # cell movement
  cell_movement!(model, itp)
end
