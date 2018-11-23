using Test
using Neurodev

function test_domain(out::Float64)
    @test (out >= -1) && (out <= 1)
end

function test_domain(out::Array{Float64})
    @test all((out .>= -1) .& (out .<= 1))
end

function test_model(m::Model)
    @test length(m.cells) < m.cfg["cells_max"]
    for cell in m.cells
        test_domain(cell.inputs)
        @test length(cell.inputs) == m.cfg["n_channels"]
        test_domain(cell.outputs)
        @test length(cell.outputs) == m.cfg["n_channels"]
        test_domain(cell.state)
        @test length(cell.state) == m.cfg["n_cell_state"]
        test_domain(cell.params)
        @test length(cell.params) == m.cfg["n_cell_params"]
        for conn in cell.conns
            test_domain(conn.inputs)
            @test length(conn.inputs) == m.cfg["n_channels"]
            test_domain(conn.outputs)
            @test length(conn.outputs) == m.cfg["n_channels"]
            test_domain(conn.state)
            @test length(conn.state) == m.cfg["n_conn_state"]
            test_domain(conn.params)
            @test length(conn.params) == m.cfg["n_conn_params"]
            @test conn.source[] === cell
            for c in eachindex(conn.cells)
                @test conn.cells[c][] in m.cells
            end
        end
    end
end

@testset "Initialization tests" begin
    cfg = Config()
    c = Controller(cfg)
    nin = rand(5:10)
    nout = rand(5:10)

    @testset "Layered" begin
        m = Model(cfg, c)
        layered_init!(m, nin, nout)
        test_model(m)
    end

    @testset "Layered single cell" begin
        m = Model(cfg, c)
        layered_init!(m, nin, nout; nhidden=1)
        test_model(m)
    end

    @testset "Random cell" begin
        m = Model(cfg, c)
        random_init!(m, nin, nout)
        test_model(m)
    end

    @testset "Random single cell" begin
        m = Model(cfg, c)
        random_init!(m, nin, nout; nhidden=1)
        test_model(m)
    end
end

@testset "Model tests" begin
    cfg = Config()
    c = Controller(cfg)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)

    @testset "Initialize" begin
        random_init!(m, nin, nout)
        test_model(m)
    end

    @testset "Input" begin
        set_input!(m, rand(nin))
        test_model(m)
    end

    @testset "Step" begin
        for i in 1:5
            step!(m)
            test_model(m)
        end
    end

    @testset "Output" begin
        outputs = get_output(m)
        test_model(m)
        test_domain(outputs)
        @test length(outputs) == nout
    end

    @testset "Reward" begin
        reward!(m, rand(nout))
        test_model(m)
    end

    @testset "Full step" begin
        for i in 1:5
            @timev outputs = step!(m, rand(nin))
            reward!(m, rand(nout))
            test_model(m)
            test_domain(outputs)
            @test length(outputs) == nout
        end
    end
end
