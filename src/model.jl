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
  ntconc::Float64
end

function Cell(id::Int64)
  pos = 1.0.+rand(length(DIMS)).*(DIMS-1.0)
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(0, id, pos, params, ctype, 0.0)
end

# type Synapse
#   c1::Cell
#   c2::Cell
#   weight::Float64
# end

type Model
  morphogens::Array{Float64}
  cells::Dict{Int64, Cell}
  # synapses::Array{Synapse}
  soma_axons::Dict{Int64, Vector{Int64}}
  synapse::Dict{Int64, Int64}
  synapse_weights::Dict{Int64, Float64}
  itp::InterpGrid
end

function Model()
  i = 1
  morphogens = zeros(Float64, [DIMS[:];N_MORPHS]...)
  param_permutes = Array{Int64}(N_MORPHS^N_PARAMS, N_PARAMS)
  i = 1
  for c = Counter([N_MORPHS for m=1:N_PARAMS])
    param_permutes[i,:] = c
    i += 1
  end
  n_cuts = ceil(Int64,N_MORPHS^(N_PARAMS/N_D))
  cuts = [linspace(1,DIMS[d],n_cuts) for d=1:N_D]
  i = 1
  cells = Dict{Int64, Cell}()
  for c = Counter([n_cuts for m=1:N_D])
    if i <= size(param_permutes)[1]
      pos = convert(Array{Float64},[cuts[d][c[d]] for d=1:N_D])
      params = vec(param_permutes[i,:])
      cells[i] = Cell(0, i, pos, params, 1, 0.0)
      i += 1
    else
      break
    end
  end
  soma_axons = Dict{Int64, Vector{Int64}}()
  synapse = Dict{Int64, Int64}()
  synapse_weights = Dict{Int64, Float64}()
  itp = InterpGrid(morphogens, BCnil, InterpQuadratic)
  Model(morphogens, cells, soma_axons, synapse, synapse_weights, itp)
end

function maxid(model::Model)
  maximum(keys(model.cells))
end

function cell_inputs(model::Model, cell::Cell)
  mid = maxid(model)
  CellInputs(cell.id/mid, cell.p_id/mid, cell.ctype, cell.params, cell.ntconc)
end

function add_cell!(model::Model, cell::Cell)
  if haskey(model.cells, cell.id)
    error("cell already in model")
  end
  model.cells[cell.id] = cell
  if cell.ctype == 3
    model.soma_axons[cell.id] = []
  elseif cell.ctype == 4
    if ~haskey(model.soma_axons, cell.p_id)
      model.soma_axons[cell.p_id] = []
    end
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
    if haskey(model.soma_axons, cell.p_id)
      soma_axons = model.soma_axons[cell.p_id]
      deleteat!(soma_axons, findin(soma_axons, cell.id))
    end
  end
  delete!(model.cells, cell.id)
end

function update_morphogens!(model::Model, cont::Controller)
  # TODO: n-dimensionality without performance hit
  for x in 1:DIMS[1]
    for y in 1:DIMS[2]
      for z in 1:DIMS[3]
        for (ckey, cell) in model.cells
          dist = evaluate(Euclidean(), cell.pos, [x y z])
          morphs = cont.morphogen_diff(dist, cell_inputs(model, cell))
          for m in 1:N_MORPHS
            model.morphogens[x,y,z,m] += morphs[m]
          end
        end
      end
    end
  end
  model.morphogens = max(0.0, model.morphogens)
  model.itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)
  nothing
  # morphogen normalization
  # sm = [size(morphs)[1:N_D]...]
  # maxmorphs = vec(maximum(model.morphogens, collect(1:N_D)))
  # for m in eachindex(maxmorphs)
  #   model.morphogens[[[1:sm[i] for i in eachindex(sm)];m]...] /= maxmorphs[m]
  # end
end

function cell_division!(model::Model, cont::Controller)
  apop = Array{Int64}(0)
  new_cells = Array{Cell}(0)
  for (ckey,cell) in model.cells
    morphogens = Array{Float64}([max(0.0,model.itp[[cell.pos[:];m]...]) for m in 1:N_MORPHS])
    if cont.division(morphogens, cell_inputs(model, cell))
      ctype = cont.child_type(morphogens, cell_inputs(model, cell))
      id = maxid(model) + length(new_cells) + 1
      params = cont.child_params(morphogens, ctype, cell_inputs(model, cell))
      pos = cell.pos + cont.child_position(morphogens, cell_inputs(model, cell)).*DIMS
      pos = min(max(pos, [1.,1.,1.]), DIMS)
      p_id = cell.id
      new_cell = Cell(p_id, id, pos, params, ctype, 0.0)
      append!(new_cells, [new_cell])
    end
    if cont.apoptosis(morphogens, cell_inputs(model, cell))
      append!(apop, [cell.id])
    end
  end
  for n in new_cells
    add_cell!(model, n)
  end
  for a in apop
    if haskey(model.cells, a) # soma deletion will cause axon deletion
      apoptosis!(model, model.cells[a])
    end
  end
  nothing
end

function synapse_update!(model::Model, cont::Controller)
  somas = filter((k,v)->v.ctype==3, model.cells)
  axons = filter((k,v)->v.ctype==4, model.cells)
  for (akey, acell) in axons
    # check for possible new synapses
    if ~haskey(model.synapse, akey)
      for (skey, scell) in somas
        if skey != acell.p_id
          dist = evaluate(Euclidean(), acell.pos, scell.pos)/mean(DIMS)
          if cont.synapse_formation(dist, cell_inputs(model, acell), cell_inputs(model, scell))
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
      smorphs = Array{Float64}([max(model.itp[[soma.pos[:];m]...],0.0) for m in 1:N_MORPHS])
      amorphs = Array{Float64}([max(model.itp[[acell.pos[:];m]...],0.0) for m in 1:N_MORPHS])
      #TODO: include reinforcement signal as reward input
      model.synapse_weights[akey] += cont.synapse_weight(0.0, smorphs, amorphs, cell_inputs(model, acell),
                                                         cell_inputs(model, soma))
    end
  end
  nothing
end

function cell_movement!(model::Model, cont::Controller)
  for (ckey,cell) in model.cells
    morphs = Vector{Float64}(N_MORPHS)
    grad = Array{Float64}(N_MORPHS, N_D)
    for m = 1:N_MORPHS
      v,g = valgrad(model.itp, [cell.pos[:];m]...)
      morphs[m] = max(v,0.0)
      grad[m,:] = g[1:N_D]
    end
    mov = cont.cell_movement(morphs, grad, cell_inputs(model, cell))
    if any(x->x>0, mov)
      mov = mov .* DIMS
      new_pos = min(max(cell.pos+mov, [1.,1.,1.]), DIMS)
      # cell.velocity = mapreduce(x->x^2, +, new_pos-cell.pos)
      cell.pos = new_pos
    end
  end
  nothing
end

function step!(model::Model, cont::Controller)
  # update morphogens
  update_morphogens!(model, cont)

  # division and apoptosis
  cell_division!(model, cont)

  # synapse formation and checking
  synapse_update!(model, cont)

  # cell movement
  cell_movement!(model, cont)
end
