using Base.Test
using NGE

@testset "Model Tests" begin

  m = Model()
  c = Controller()

  @testset "Step" begin
    step!(m, c)
    @test 0 == 0 # just test for errors
  end

  @testset "Add cells" begin
    ncells = length(m.cells)
    g = Cell()
    add_cell!(m, g)
    @test length(findin(m.cells, [g])) == 1
    @test length(m.cells) == ncells + 1
  end

  @testset "Apoptosis" begin
    ncells = length(m.cells)
    c1 = m.cells[1]
    apoptosis!(m, c1)
    @test length(findin(m.cells, [c1])) == 0
    @test length(m.cells) == ncells - 1
    c1 = Cell()
    c2 = Cell()
    add_cell!(m, c1)
    add_cell!(m, c2)
    nsyns = length(m.synapses)
    add_synapse!(m, c1, c2)
    @test length(m.synapses) == nsyns + 1
    @test length(find(x->x.c1 == c1 && x.c2 == c2, m.synapses)) == 1
    apoptosis!(m, c1)
    @test length(findin(m.cells, [c1])) == 0
    @test length(m.synapses) == nsyns
  end

  @testset "Update morphogens" begin
    c.morphogen_diff = (morphogens::Vector{Float64}, dist::Vector{Float64}, cell::CellInputs) -> randn(N_MORPHS)
    update_morphogens!(m, c)
    @test any(m.morphogens .> 0)
    @test ~any(m.morphogens .< 0)
  end

  @testset "Cell division" begin
    c.division = (morphogens::Vector{Float64}, cell::CellInputs) -> true
    add_cell!(m, Cell())
    cells = deepcopy(m.cells)
    cell_division!(m, c)
    @test cells != m.cells
    @test length(m.cells) > length(cells)
  end

  @testset "Cell movement" begin
    c.cell_movement = (morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs) -> rand(N_D)
    pos = map(x->x.pos, m.cells)
    cell_movement!(m, c)
    newpos = map(x->x.pos, m.cells)
    @test newpos != pos
    @test all(x->all(x.>=1.), newpos)
    @test all(x->all(x.<=DIMS), newpos)
  end

  @testset "Synapse formation" begin
    nsyns = length(m.synapses)
    c1 = Cell()
    c2 = Cell()
    add_cell!(m, c1)
    add_cell!(m, c2)

    c.synapse_formation = (dist::Float64, acell::CellInputs, bcell::CellInputs) -> true
    c.synapse_weight =
      (a_morphs::Vector{Float64}, b_morphs::Vector{Float64}, acell::CellInputs, bcell::CellInputs) -> 0.0

    # form synapses
    synapse_update!(m, c)

    # neural cells should be connected with new neural weights
    @test length(m.synapses) > nsyns
    @test length(find(x->x.c1 == c1 && x.c2 == c2, m.synapses)) == 1
    @test length(find(x->x.c1 == c2 && x.c2 == c1, m.synapses)) == 1
    @test all(x->m.synapses[x].weight == 0.0,
              find(x->x.c1 == c1 || x.c1 == c2 || x.c2 == c1 || x.c2 == c2, m.synapses))

    c.synapse_weight =
      (a_morphs::Vector{Float64}, b_morphs::Vector{Float64}, acell::CellInputs, bcell::CellInputs) -> randn()

    # cache the current weights and update
    weights = deepcopy(map(s->s.weight, m.synapses))
    nsyns = length(m.synapses)
    synapse_update!(m, c)

    @test length(m.synapses) == nsyns
    @test weights != map(s->s.weight, m.synapses)
  end

  @testset "Firing" begin
    # ensure enough cells
    while length(m.cells) < 10
      add_cell!(m, Cell())
    end

    # ensure synaptic connections with nonzero weights
    c.synapse_formation = (dist::Float64, acell::CellInputs, bcell::CellInputs) -> true
    synapse_update!(m, c)
    for s in filter(x->x.weight == 0.0, m.synapses)
      s.weight = randn()
    end

    ntconcs = deepcopy(map(x->x.ntconc, m.cells))

    c.nt_output = (input::Float64, cell::CellInputs) -> randn()
    c.nt_update = (input::Float64, output::Float64, cell::CellInputs) -> input

    fire!(m, c) # fire into next timestep inputs
    fire!(m, c) # nt update each cell

    @test ntconcs != map(x->x.ntconc, m.cells)
    @test any(x->x.ntin != 0.0, m.cells)

    # set to null output and compare concentration
    ntins = deepcopy(map(x->x.ntin, m.cells))
    ntconcs = deepcopy(map(x->x.ntconc, m.cells))
    c.nt_output = (input::Float64, cell::CellInputs) -> 0.0

    fire!(m, c)

    @test ntins + ntconcs == map(x->x.ntconc, m.cells)
    @test all(x->x.ntin == 0.0, m.cells)

    # multiple firing doesn't change, as nt_out is 0.0
    fire!(m, c)

    @test ntins + ntconcs == map(x->x.ntconc, m.cells)
    @test all(x->x.ntin == 0.0, m.cells)
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
