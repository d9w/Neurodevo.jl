using LightGraphs
using Grid

include("controller.jl")
include("constants.jl")

type Cell
  pos::Vector{Float64}
  params::Vector{Int64}
  ctype::Int64
  ntconc::Float64
  ntin::Float64
end

function Cell()
  pos = 1.0.+rand(length(DIMS)).*(DIMS-1.0)
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(pos, params, ctype, 0.0, 0.0)
end

type Synapse
  c1::Cell # reference
  c2::Cell
  weight::Float64
end

type Model
  morphogens::Array{Float64}
  cells::Dict{Int64, Cell}
  synapses::Array{Synapse}
  itp::InterpGrid
end

function Model()
  i = 1
  morphogens = zeros(Float64, [DIMS[:];N_MORPHS]...)
  # TODO: this is still making 1 of each param type
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
      cells[i] = Cell(i, pos, params, 1, 0.0, 0.0)
      i += 1
    else
      break
    end
  end
  synapses= Array{Synapse}(0)
  itp = InterpGrid(morphogens, BCnil, InterpQuadratic)
  Model(morphogens, cells, synapses, itp)
end

function posm(cell_pos::Vector{Float64}, m::Int64)
  # TODO: n-dimensionality without performance hit
  [cell_pos[1], cell_pos[2], cell_pos[3], m]
end

function CellInputs(cell::Cell)
  CellInputs(cell.ctype, cell.params, cell.ntconc)
end

function add_cell!(model::Model, cell::Cell)
  if findfirst(model.cells, cell)
    error("cell already in model")
  end
  append!(model.cells, [cell])
end

function apoptosis!(model::Model, cell::Cell)
  filter!(s->((s.c1 != cell) && (s.c2 != cell)), model.synapses)
  deleteat!(model.cells, findfirst(model.cells, cell))
end

function update_morphogens!(model::Model, cont::Controller)
  # TODO: n-dimensionality without performance hit
  for x in 1:DIMS[1]
    for y in 1:DIMS[2]
      for z in 1:DIMS[3]
        for cell in model.cells
          dvec = [x y z] - cell.pos
          # TODO: reward
          morphs = cont.morphogen_diff(model.morphogens[x,y,z,m], dvec, CellInputs(cell), 0.0)
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
  apop_cells = Array{Cell}(0)
  new_cells = Array{Cell}(0)
  for cell in model.cells
    # morphogens = Array{Float64}([max(0.0,model.itp[[cell.pos[:];m]...]) for m in 1:N_MORPHS])
    morphogens = zeros(N_MORPHS) # until using models with cell division
    if cont.division(morphogens, CellInputs(cell))
      ctype = cont.child_type(morphogens, CellInputs(cell))
      params = cont.child_params(morphogens, ctype, CellInputs(cell))
      pos = cell.pos + cont.child_position(morphogens, CellInputs(cell)).*DIMS
      pos = min(max(pos, [1.,1.,1.]), DIMS)
      new_cell = Cell(pos, params, ctype, 0.0, 0.0)
      push!(new_cells, new_cell)
    end
    if cont.apoptosis(morphogens, CellInputs(cell))
      push!(apop_cells, cell)
    end
  end
  for cell in new_cells
    add_cell!(model, cell)
  end
  for cell in apop_cells
    apoptosis!(model, cell)
  end
  nothing
end

function synapse_update!(model::Model, cont::Controller)
  bsynapses = zeros(Bool, length(model.cells), length(model.cells))
  for s in model.synapses
    # TODO: eval performance
    bsynapses[findfirst(model.cells, s.c1), findfirst(model.cells, s.c2)] = true
  end

  for c1 in eachindex(model.cells)
    for c2 in eachindex(model.cells)
      if c1 != c2 && ~former_synapse[c1, c2]
        c1cell = model.cells[c1]
        c2cell = model.cells[c2]
        dist = euclidean(c1cell.pos, c2cell.pos)/mean(DIMS)
        if cont.synapse_formation(dist, CellInputs(c1cell), CellInputs(c2cell))
          append!(model.synapses, [Synapse(c1cell, c2cell, 0.0)])
        end
      end
    end
  end

  for s in model.synapses
    a_morphs = Array{Float64}([max(model.itp[posm(s.c1.pos,m)...],0.0) for m in 1:N_MORPHS])
    b_morphs = Array{Float64}([max(model.itp[posm(s.c2.pos,m)...],0.0) for m in 1:N_MORPHS])
    s.weight += cont.synapse_weight(a_morphs, b_morphs, CellInputs(s.c1), CellInputs(s.c2))
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
    mov = cont.cell_movement(morphs, grad, CellInputs(cell))
    if any(x->x>0, mov)
      mov = mov .* DIMS
      new_pos = min(max(cell.pos+mov, [1.,1.,1.]), DIMS)
      # cell.velocity = mapreduce(x->x^2, +, new_pos-cell.pos)
      cell.pos = new_pos
    end
  end
  nothing
end

function fire!(model::Model, cont::Controller)
  outputs = zeroes(length(model.cells))
  for c in eachindex(model.cells)
    cell = model.cells[c]
    outputs[c] = cont.nt_output(cell.nt_in, CellInputs(cell))
    cell.ntconc += cont.nt_update(cell.nt_in, outputs[c], CellInputs(cell))
    cell.nt_in = 0.0
  end

  for c in eachindex(model.cells)
    if outputs[c] != 0.0
      for s in filter(s->s.c1==c, model.synapses)
        model.cells[s.c2].nt_in += outputs[c] * s.weight
      end
    end
  end
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
