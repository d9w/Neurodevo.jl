using Base.Test
using NGE

@testset "Model Tests" begin

  m = Model()
  c = Controller()

  @testset "Add cells" begin
    ncells = length(m.cells)
    g = Cell()
    add_cell!(m, g)
    @test length(findin(m.cells, [g])) == 1
    @test length(m.cells) == ncells + 1
  end

  @testset "Apoptosis" begin
    c = m.cells[1]
    apoptosis!(m, c)
    @test length(findin(m.cells, [c])) == 0
    @test length(m.cells) == 0
    c1 = Cell()
    c2 = Cell()
    add_cell!(m, c1)
    add_cell!(m, c2)
    @test length(m.synapses) == 0
    add_synapse!(m, c1, c2)
    @test length(m.synapses) == 1
    @test m.synapses[1].c1 == c1
    @test m.synapses[1].c2 == c2
    apoptosis!(m, c1)
    @test length(findin(m.cells, [c1])) == 0
    @test length(m.synapses) == 0
  end

  @testset "Update morphogens" begin
    cont = Controller()
    update_morphogens!(m, cont)
    @test any(m.morphogens .> 0)
    @test ~any(m.morphogens .< 0)
  end

  @testset "Cell division" begin
    cont.division = (morphogens::Vector{Float64}, cell::CellInputs) -> true
    cells = deepcopy(m.cells)
    cell_division!(m, cont)
    @test cells != m.cells
    @test length(m.cells) > length(cells)
  end

  @testset "Cell movement" begin
    pos = map(x->x.pos, m.cells)
    update_morphogens!(m, cont)
    cell_movement!(m, cont)
    newpos = map(x->x.pos, m.cells)
    @test newpos != pos
    @test all(x->all(x.>=1.), newpos)
    @test all(x->all(x.<=DIMS), newpos)
  end

  @testset "Synapse formation" begin
    c1 = Cell()
    c2 = Cell()
    add_cell!(m, c1)
    add_cell!(m, c2)

    # force synapse formation and survival
    cont = Controller()
    cont.synapse_formation = (dist::Float64, acell::CellInputs, bcell::CellInputs) -> true

    # form synapses
    synapse_update!(m, cont)

    # neural cells should be connected with new neural weights
    @test length(m.synapses) == 2
    @test m.synapses[1].c1 == c1 || m.synapses[1].c2 == c1
    @test m.synapses[1].c1 == c2 || m.synapses[1].c2 == c2
    @test m.synapses[2].c1 == c1 || m.synapses[2].c2 == c1
    @test m.synapses[2].c1 == c2 || m.synapses[2].c2 == c2

    weights = deepcopy(map(s->s.weight, m.synapses))

    # cache the current weights and update
    synapse_update!(m, cont)
    @test length(m.synapses) == 2
    @test weights != map(s->s.weight, m.synapses)
  end

end

# function profile_model()
#   model = Model()
#   cont = Controller()
#   sumt = 0.0
#   suma = 0
#   n_steps = 10
#   for c = 1:n_steps
#     tic()
#     suma += @allocated step!(model, cont)
#     sumt += toq()
#     println([sumt/n_steps, suma/n_steps, length(model.cells), length(model.synapses)])
#   end
# end
