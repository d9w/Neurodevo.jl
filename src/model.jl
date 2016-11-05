using Grid
using Distances

function Cell()
  pos = 1.0.+rand(length(DIMS)).*(DIMS-1.0)
  params = rand(1:N_MORPHS, N_PARAMS)
  ctype = 1 # neural stem cell
  Cell(pos, params, ctype, 0.0, 0.0)
end

function posm(cell_pos::Vector{Float64}, m::Int64)
  # TODO: n-dimensionality without performance hit
  [cell_pos[1], cell_pos[2], cell_pos[3], m]
end

function CellInputs(cell::Cell)
  CellInputs(cell.ctype, cell.params, cell.ntconc)
end

function add_cell!(model::Model, cell::Cell)
  if findfirst(model.cells, cell) > 0
    error("cell already in model")
  end
  push!(model.cells, cell)
end

function add_synapse!(model::Model, c1::Cell, c2::Cell)
  for c in [c1, c2]
    if findfirst(model.cells, c) == 0
      add_cell!(model, c)
    end
  end
  s = Synapse(c1, c2, 0.0)
  push!(model.synapses, s)
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
          dvec = [x, y, z] - cell.pos
          # TODO: reward
          morphs = cont.morphogen_diff([model.morphogens[x,y,z,m] for m in 1:N_MORPHS], dvec, CellInputs(cell))
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
      if c1 != c2 && ~bsynapses[c1, c2]
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
  for cell in model.cells
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
