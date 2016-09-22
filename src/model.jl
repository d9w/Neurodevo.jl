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
  velocity::Float64
end

function Cell(id::Int64)
  pos = 1.0.+rand(length(DIMS)).*(DIMS-1.0)
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(0, id, pos, params, ctype, 0.0)
end

type Model
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
  Model(morphogens, grid, cells, soma_axons, synapse, synapse_weights, synapse_graph)
end

function maxid(model::Model)
  maximum(keys(model.cells))
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
    if haskey(model.soma_axons, cell.p_id)
      soma_axons = model.soma_axons[cell.p_id]
      deleteat!(soma_axons, findin(soma_axons, cell.id))
    end
  end
  delete!(model.cells, cell.id)
end

function update_morphogens!(model::Model, cont::Controller)
  # model.morphogens *= MORPHOGEN_DECAY
  for i in eachindex(model.grid[:,1])
    pos = vec(model.grid[i,:])
    for (ckey, cell) in model.cells
      for m in 1:N_MORPHS
        morph = model.morphogens[[pos;m]...]
        dist = evaluate(Euclidean(), cell.pos, pos)
        morph += cont.morphogen_diff(N_MORPHS, m, cell.ctype, cell.params, dist)
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
  apop = Array{Int64}(0)
  new_cells = Array{Cell}(0)
  for (ckey,cell) in model.cells
    morphogens = Array{Float64}([itp[[cell.pos[:];m]...] for m in 1:N_MORPHS])
    branch = cont.division(morphogens, cell.ctype, cell.params, cell.velocity)
    if branch > 0
      branch_dir = branch == 1
      ctype = cont.child_type(morphogens, cell.ctype, cell.params, branch_dir)
      if ctype == 0
        append!(apop, [cell.id])
      else
        id = maxid(model) + length(new_cells) + 1
        params = cont.child_params(morphogens, cell.ctype, cell.params, ctype)
        pos = cell.pos + cont.child_position(morphogens, cell.ctype, cell.params).*DIMS
        pos = min(max(pos, [1.,1.,1.]), DIMS)
        p_id = cell.id
        if cell.ctype == 4
          p_id = cell.p_id
        end
        new_cell = Cell(p_id, id, pos, params, ctype, 0.0)
        append!(new_cells, [new_cell])
      end
    end
  end
  for a in apop
    if haskey(model.cells, a) # soma deletion will cause axon deletion
      apoptosis!(model, model.cells[a])
    end
  end
  for n in new_cells
    add_cell!(model, n)
  end
end

function synapse_update!(model::Model, itp::Grid.InterpGrid, cont::Controller)
  somas = filter((k,v)->v.ctype==3, model.cells)
  axons = filter((k,v)->v.ctype==4, model.cells)
  for (akey, acell) in axons
    # check for possible new synapses
    if ~haskey(model.synapse, akey)
      for (skey, scell) in somas
        if skey != acell.p_id
          dist = evaluate(Euclidean(), acell.pos, scell.pos)
          if cont.synapse_formation(dist)
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
      #TODO: include reinforcement signal as reward input
      model.synapse_weights[akey] += cont.synapse_weight(smorphs, amorphs, soma.params, acell.params, 0.0)
      if ~cont.synapse_survival(model.synapse_weights[akey])
        apoptosis!(model, acell)
      end
    end
  end
end

function cell_movement!(model::Model, itp::Grid.InterpGrid, cont::Controller)
  for (ckey,cell) in model.cells
    morphs = Vector{Float64}(N_MORPHS)
    grad = Array{Float64}(N_MORPHS, N_D)
    for m = 1:N_MORPHS
      v,g = valgrad(itp, [cell.pos[:];m]...)
      morphs[m] = v
      grad[m,:] = g[1:N_D]
    end
    if ~in(ckey, [collect(keys(model.synapse));collect(values(model.synapse))])
      mov = cont.cell_movement(morphs, grad, cell.ctype, cell.params)
      if any(x->x>0, mov)
        mov = mov .* DIMS
        cell.velocity = mapreduce(x->x^2, +, mov)
        cell.pos = min(max(cell.pos+mov, [1.,1.,1.]), DIMS)
      end
    end
  end
end

function step!(model::Model, cont::Controller)
  # update morphogens
  update_morphogens!(model, cont)
  itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)

  # division and apoptosis
  cell_division!(model, itp, cont)

  # synapse formation and checking
  synapse_update!(model, itp, cont)

  # cell movement
  cell_movement!(model, itp, cont)
end
