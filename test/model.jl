using StatsBase
using Gadfly
using Colors

include("../src/model.jl")

function add_cells()
  m = Model()
  g = Cell(maxid(m)+1)
  g.ctype = 2
  add_cell!(m, g)
  @assert haskey(m.cells, g.id)
  @assert maxid(m)==g.id
  soma = Cell(maxid(m)+1)
  soma.p_id = 1
  soma.ctype = 3
  add_cell!(m, soma)
  @assert haskey(m.soma_axons, soma.id)
  axon = Cell(maxid(m)+1)
  axon.p_id = soma.id
  axon.ctype = 4
  add_cell!(m, axon)
  @assert haskey(m.synapse, axon.id)
  @assert haskey(m.synapse_weight, axon.id)
  @assert axon.id in m.soma_axons[soma.id]
  m
end

function apoptosis()
  m = add_cells()
  apoptosis!(m, m.cells[1])
  @assert ~haskey(m.cells,1)
  soma = first(keys(filter((k,v)->v.ctype==3, m.cells)))
  axon = first(keys(filter((k,v)->v.ctype==4, m.cells)))
  apoptosis!(m, m.cells[axon])
  @assert ~haskey(m.cells, axon)
  @assert ~haskey(m.synapse, axon)
  @assert ~haskey(m.synapse_weight, axon)
  @assert ~(axon in m.soma_axons[soma])
  @assert mapreduce(x->axon in x[2], +, m.soma_axons) == 0
  new_axon = Cell(maxid(m)+1)
  new_axon.p_id = soma
  new_axon.ctype = 4
  add_cell!(m, new_axon)
  apoptosis!(m, m.cells[soma])
  @assert ~haskey(m.cells, soma)
  @assert ~haskey(m.soma_axons, soma)
  @assert ~haskey(m.cells, new_axon.id)
  @assert ~haskey(m.synapse, new_axon.id)
  @assert ~haskey(m.synapse_weight, new_axon.id)
  @assert mapreduce(x->new_axon.id in x[2], +, m.soma_axons) == 0
end

function update_morphogens()
  m = Model()
  cont = Controller()
  update_morphogens!(m, cont)
  @assert any(m.morphogens .> 0)
  @assert ~any(m.morphogens .< 0)
end

function cell_division()
  m = Model()
  cont = Controller()
  # force division decision
  cont.division = (a, b, c, d) -> true
  itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)
  cells = deepcopy(collect(values(m.cells)))
  cell_division!(m, itp, cont)
  @assert cells != collect(values(m.cells))
  newcells = filter((k,v)->v.p_id!=0,m.cells)
  @assert length(newcells) > 0
end

function cell_movement()
  m = Model()
  cont = Controller()
  pos = map(x->x[2].pos, m.cells)
  update_morphogens!(m, cont)
  itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)
  cell_movement!(model, itp, cont)
  newpos = map(x->x[2].pos, m.cells)
  @assert newpos != pos
  @assert all(x->all(x.>0), newpos)
  @assert all(x->all(x.<DIMS), newpos)
end

function synapse_update()
  m = Model()
  axons = Array{Int64}(2)
  somas = Array{Int64}(2)
  for i=1:2
    # add two neural cells
    soma = Cell(maxid(m)+1)
    soma.ctype = 3
    add_cell!(m, soma)
    somas[i] = soma.id
    axon = Cell(maxid(m)+1)
    axon.p_id = soma.id
    axon.ctype = 4
    add_cell!(m, axon)
    axons[i] = axon.id
  end

  # force synapse formation and survival
  cont = Controller()
  cont.synapse_formation = (a,b) -> true
  cont.synapse_survival = (a,b) -> true

  # cache the current weights and update
  weights = deepcopy(values(m.synapse_weights))
  itp = InterpGrid(model.morphogens, BCnil, InterpQuadratic)
  synapse_update!(m, itp, cont)

  # neural cells should be connected with new neural weights
  @assert m.synapse[axons[1]] == somas[2]
  @assert m.synapse[axons[2]] == somas[1]
  @assert weights != values(m.synapse_weights)
end

function handwritten_rules()
  model = Model()
  for i = 1:10
    step!(model)
  end
end
