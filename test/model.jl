using Base.Test
using StatsBase

include("../src/model.jl")
include("../src/graph.jl")

function add_cells()
  m = Model()
  g = Cell(maxid(m)+1)
  g.ctype = 2
  add_cell!(m, g)
  @test haskey(m.cells, g.id)
  @test maxid(m)==g.id
  soma = Cell(maxid(m)+1)
  soma.p_id = 1
  soma.ctype = 3
  add_cell!(m, soma)
  @test haskey(m.soma_axons, soma.id)
  axon = Cell(maxid(m)+1)
  axon.p_id = soma.id
  axon.ctype = 4
  add_cell!(m, axon)
  @test axon.id in m.soma_axons[soma.id]
  m
end

function apoptosis()
  m = add_cells()
  apoptosis!(m, m.cells[1])
  @test ~haskey(m.cells,1)
  soma = first(keys(filter((k,v)->v.ctype==3, m.cells)))
  axon = first(keys(filter((k,v)->v.ctype==4, m.cells)))
  apoptosis!(m, m.cells[axon])
  @test ~haskey(m.cells, axon)
  @test ~haskey(m.synapse, axon)
  @test ~haskey(m.synapse_weights, axon)
  @test ~(axon in m.soma_axons[soma])
  @test mapreduce(x->axon in x[2], +, m.soma_axons) == 0
  new_axon = Cell(maxid(m)+1)
  new_axon.p_id = soma
  new_axon.ctype = 4
  add_cell!(m, new_axon)
  apoptosis!(m, m.cells[soma])
  @test ~haskey(m.cells, soma)
  @test ~haskey(m.soma_axons, soma)
  @test ~haskey(m.cells, new_axon.id)
  @test ~haskey(m.synapse, new_axon.id)
  @test ~haskey(m.synapse_weights, new_axon.id)
  @test isempty(m.soma_axons) || (mapreduce(x->new_axon.id in x[2], +, m.soma_axons) == 0)
end

function update_morphogens()
  m = Model()
  cont = Controller()
  update_morphogens!(m, cont)
  @test any(m.morphogens .> 0)
  @test ~any(m.morphogens .< 0)
end

function cell_division()
  m = Model()
  cont = Controller()
  # force division decision
  cont.division = (morphogens::Vector{Float64}, cell::CellInputs) -> true
  cells = deepcopy(collect(values(m.cells)))
  cell_division!(m, cont)
  @test cells != collect(values(m.cells))
  newcells = filter((k,v)->v.p_id!=0,m.cells)
  @test length(newcells) > 0
end

function cell_movement()
  m = Model()
  cont = Controller()
  pos = map(x->x[2].pos, m.cells)
  update_morphogens!(m, cont)
  cell_movement!(m, cont)
  newpos = map(x->x[2].pos, m.cells)
  @test newpos != pos
  @test all(x->all(x.>=1.), newpos)
  @test all(x->all(x.<=DIMS), newpos)
end

function synapse_update()
  m = Model()
  for (ckey, cell) in m.cells
    delete!(m.cells, ckey)
  end
  axons = Array{Int64}(2)
  somas = Array{Int64}(2)
  id = 1
  for i=1:2
    # add two neural cells
    soma = Cell(id)
    soma.ctype = 3
    somas[i] = soma.id
    add_cell!(m, soma)
    id+=1
    axon = Cell(id)
    axon.p_id = soma.id
    axon.ctype = 4
    axons[i] = axon.id
    add_cell!(m, axon)
    id+=1
  end

  # force synapse formation and survival
  cont = Controller()
  cont.synapse_formation = (dist::Float64, acell::CellInputs, bcell::CellInputs) -> true

  # cache the current weights and update
  weights = deepcopy(values(m.synapse_weights))
  synapse_update!(m, cont)

  # neural cells should be connected with new neural weights
  @test m.synapse[axons[1]] == somas[2]
  @test m.synapse[axons[2]] == somas[1]
  @test weights != values(m.synapse_weights)
end

function test_all()
  print("add_cells")
  add_cells();
  print_with_color(:green, "...[passed]\n")
  print("apoptosis")
  apoptosis()
  print_with_color(:green, "...[passed]\n")
  print("update_morphogens")
  update_morphogens()
  print_with_color(:green, "...[passed]\n")
  print("cell_division")
  cell_division()
  print_with_color(:green, "...[passed]\n")
  print("cell_movement")
  cell_movement()
  print_with_color(:green, "...[passed]\n")
  print("synapse_update")
  synapse_update()
  print_with_color(:green, "...[passed]\n")
end

function profile_model()
  model = Model()
  cont = random_controller()
  sumt = 0.0
  suma = 0
  n_steps = 10
  for c = 1:n_steps
    tic()
    suma += @allocated step!(model, cont)
    sumt += toq()
    println([sumt,suma]/n_steps)
  end
end
