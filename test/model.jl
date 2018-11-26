using Test
using Neurodev

function test_model(m::Model, nin::Int64, nout::Int64, nhidden::Int64)
    @test length(m.cells) < m.cfg["cells_max"]
    for cell in m.cells
        @test all(cell.inputs .>= -1.0)
        @test all(cell.inputs .<= 1.0)
        @test length(cell.inputs) == m.cfg["n_channels"]
        @test all(cell.outputs .>= -1.0)
        @test all(cell.outputs .<= 1.0)
        @test length(cell.outputs) == m.cfg["n_channels"]
        @test all(cell.state .>= -1.0)
        @test all(cell.state .<= 1.0)
        @test length(cell.state) == m.cfg["n_cell_state"]
        @test all(cell.params .>= -1.0)
        @test all(cell.params .<= 1.0)
        @test length(cell.params) == m.cfg["n_cell_params"]
        for conn in cell.conns
            @test all(conn.state .>= -1.0)
            @test all(conn.state .<= 1.0)
            @test length(conn.state) == m.cfg["n_conn_state"]
            @test all(conn.params .>= -1.0)
            @test all(conn.params .<= 1.0)
            @test length(conn.params) == m.cfg["n_conn_params"]
            @test conn.source[] === cell
            @test conn.dest[] in m.cells
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
        test_model(m, nin, nout, nin)
    end

    @testset "Layered single cell" begin
        m = Model(cfg, c)
        layered_init!(m, nin, nout; nhidden=1)
        test_model(m, nin, nout, 1)
    end

    @testset "Random cell" begin
        m = Model(cfg, c)
        random_init!(m, nin, nout)
        test_model(m, nin, nout, nin)
    end

    @testset "Random single cell" begin
        m = Model(cfg, c)
        random_init!(m, nin, nout; nhidden=1)
        test_model(m, nin, nout, 1)
    end
end

@testset "Model tests" begin
    cfg = Config()
    c = Controller(cfg)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)
    nhidden = nin

    @testset "Initialize" begin
        random_init!(m, nin, nout)
        test_model(m, nin, nout, nin)
        @test length(m.cells) == nin + nout + nhidden
        @test length(m.inputs) == nin
        @test length(m.outputs) == nout
    end

    @testset "Input" begin
        set_input!(m, rand(nin))
        test_model(m, nin, nout, nin)
    end

    @testset "Step" begin
        for i in 1:5
            step!(m)
            test_model(m, nin, nout, nin)
        end
    end

    @testset "Output" begin
        outputs = get_output(m)
        test_model(m, nin, nout, nin)
        @test all(outputs .>= -1.0)
        @test all(outputs .<= 1.0)
        @test length(outputs) == nout
    end

    @testset "Reward" begin
        reward!(m, rand(nout))
        test_model(m, nin, nout, nin)
    end

    @testset "Full step" begin
        for i in 1:5
            println("Base controller step $i")
            @time outputs = step!(m, rand(nin))
            reward!(m, rand(nout))
            test_model(m, nin, nout, nin)
            @test all(outputs .>= -1.0)
            @test all(outputs .<= 1.0)
            @test length(outputs) == nout
        end
    end
end

function controller_test(cont_constructor::Function; cfg::Dict=Config())
    c = cont_constructor(cfg)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)
    random_init!(m, nin, nout)

    @testset "Full step" begin
        for i in 1:5
            println(string(cont_constructor, " step ", i))
            @time outputs = step!(m, rand(nin))
            reward!(m, rand(nout))
            test_model(m, nin, nout, nin)
            @test all(outputs .>= -1.0)
            @test all(outputs .<= 1.0)
            @test length(outputs) == nout
        end
    end
end

@testset "Random controller" begin
    controller_test(rand_controller)
end

@testset "Static controller" begin
    controller_test(static_controller;
                    cfg = Config(Config(), "cfg/static.yaml"))
end

@testset "SNN controller" begin
    controller_test(snn_controller;
                    cfg = Config(Config(), "cfg/snn.yaml"))
end
