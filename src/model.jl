using LightGraphs
using Distances
using Grid

include("controller.jl")
include("constants.jl")

type Cell
  p_id::Int64
  id::Int64
  pos::Vector{Float64}
  params::Vector{Int64}
  ctype::Int64
end

function Cell(id::Int64)
  pos = 1+rand(length(DIMS)).*(DIMS-1.0)
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(0, id, pos, params, ctype)
end

type Model
  morphogens::Array{Float64}
  grid::Array{Int64}
  cells::Dict{Int64, Cell}
  soma_axons::Dict{Int64, Vector{Int64}}
  synapse::Dict{Int64, Int64}
  synapse_weights::Dict{Int64, Float64}
  synapse_graph::DiGraph
  entropy::Int64
end

function Model()
  i = 1
  grid = Array{Int64}(*(DIMS...),length(DIMS))
  for c = Counter(DIMS)
    grid[i,:] = c
    i += 1
  end
  morphogens = zeros(Float64, [DIMS[:];N_MORPHS]...)
  cells = Dict{Int64, Cell}()
  for i=1:N_NSC
    cell = Cell(i)
    cells[i] = cell
  end
  soma_axons = Dict{Int64, Vector{Int64}}()
  synapse = Dict{Int64, Int64}()
  synapse_weights = Dict{Int64, Float64}()
  synapse_graph = DiGraph()
  entropy = 0
  Model(morphogens, grid, cells, soma_axons, synapse, synapse_weights, synapse_graph, entropy)
end

function maxid(model::Model)
  maximum(keys(model.cells))
end

function make_bias!(model::Model)
  m.entropy += 1
  rand()
end

function add_cell!(model::Model, cell::Cell)
  if haskey(model.cells, cell.id)
    error("cell already in model")
  end
  model.cells[cell.id] = cell
  if cell.ctype == 3
    model.soma_axons[cell.id] = []
  elseif cell.ctype == 4
    push!(model.soma_axons[cell.p_id], cell.id)
  end
end

function apoptosis!(model::Model, cell::Cell)
  if cell.ctype == 3
    for a in model.soma_axons[cell.id]
      delete!(model.cells, a)
      delete!(model.synapse, a)
      delete!(model.synapse_weights, a)
    end
    delete!(model.soma_axons, cell.id)
    for s in keys(model.synapse)
      if model.synapse[s] == cell.id
        delete!(model.synapse, s)
        delete!(model.synapse_weights, s)
      end
    end
  elseif cell.ctype == 4
    delete!(model.synapse, cell.id)
    delete!(model.synapse_weights, cell.id)
    soma_id = cell.p_id
    if cell.p_id in keys(model.soma_axons)
      soma_axons = model.soma_axons[cell.p_id]
      deleteat!(soma_axons, findin(soma_axons, cell.id))
    end
  end
  delete!(model.cells, cell.id)
end

function update_morphogens!(model::Model, cont::Controller)
  # model.morphogens *= MORPHOGEN_DECAY
  bias = make_bias!(model)
  for i in eachindex(model.grid[:,1])
    pos = vec(model.grid[i,:])
    for (ckey, cell) in model.cells
      for m in 1:N_MORPHS
        morph = model.morphogens[[pos;m]...]
        dist = evaluate(Euclidean(), cell.pos, pos)
        morph += cont.morphogen_diff(N_MORPHS, m, cell.ctype, cell.params, dist, bias)
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

function cell_division!(model::Model, itp::Grid.InterpGrid, cont::Controller)
  bias = make_bias!(model)
  for (ckey,cell) in model.cells
    morphogens = Array{Float64}([itp[[cell.pos[:];m]...] for m in 1:N_MORPHS])
    if cont.division(morphogens, cell.ctype, cell.params, bias)
      branch_dir = cont.child_branch(morphogens, cell.ctype, cell.params, bias)
      ctype = cont.child_type(morphogens, cell.ctype, cell.params, branch_dir, bias)
      if ctype == 0
        apoptosis!(model, cell)
      else
        id = maxid(model) + 1
        params = cont.child_params(morphogens, cell.ctype, cell.params, ctype, bias)
        pos = cell.pos + cont.child_position(morphogens, cell.ctype, cell.params, convert(Array{Float64},DIMS), bias)
        pos = min(max(pos, [0.,0.,0.]), DIMS)
        p_id = cell.id
        if cell.ctype == 4
          p_id = cell.p_id
        end
        new_cell = Cell(p_id, id, pos, params, ctype)
        add_cell!(model, new_cell)
      end
    end
  end
end

function synapse_update!(model::Model, itp::Grid.InterpGrid, cont::Controller)
  bias = make_bias!(model)
  somas = filter((k,v)->v.ctype==3, model.cells)
  axons = filter((k,v)->v.ctype==4, model.cells)
  for (akey, acell) in axons
    # check for possible new synapses
    if ~haskey(model.synapse, akey)
      for (skey, scell) in somas
        if skey != acell.p_id
          dist = evaluate(Euclidean(), acell.pos, scell.pos)
          if cont.synapse_formation(dist, bias)
            model.synapse[akey] = skey
            model.synapse_weights[akey] = 0.0
            break
          end
        end
      end
    end
    # potentially with new synapse, update weight
    if haskey(model.synapse, akey)
      soma = model.cells[model.synapse[akey]]
      smorphs = convert(Array{Float64},[itp[[soma.pos[:];m]...] for m in 1:N_MORPHS])
      amorphs = convert(Array{Float64},[itp[[acell.pos[:];m]...] for m in 1:N_MORPHS])
      model.synapse_weights[akey] += cont.synapse_weight(smorphs, amorphs, soma.params, acell.params, bias)
      if ~cont.synapse_survival(model.synapse_weights[akey], bias)
        apoptosis!(model, acell)
      end
    end
  end
end

function cell_movement!(model::Model, itp::Grid.InterpGrid, cont::Controller)
  bias = make_bias!(model)
  for (ckey,cell) in model.cells
    morphs = Vector{Float64}(N_MORPHS)
    grad = Array{Float64}(N_MORPHS, N_D)
    for m = 1:N_MORPHS
      v,g = valgrad(itp, [cell.pos[:];m]...)
      morphs[m] = v
      grad[m,:] = g[1:N_D]
    end
    if ~in(ckey, [collect(keys(model.synapse));collect(values(model.synapse))])
      cell.pos += cont.cell_movement(morphs, grad, cell.ctype, cell.params, convert(Array{Float64}, DIMS), bias)
      cell.pos = min(max(cell.pos, [0.,0.,0.]), DIMS)
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
