function test_model(m::Model)
    @test length(m.cells) <= m.cfg["cells_max"]
    nconns = 0
    for cell in m.cells
        nconns += length(cell.conns)
        @test cell.input >= -1.0
        @test cell.input <= 1.0
        @test cell.output >= -1.0
        @test cell.output <= 1.0
        @test all(cell.state .>= -1.0)
        @test all(cell.state .<= 1.0)
        @test length(cell.state) == m.cfg["n_cell_state"]
        @test all(cell.params .>= -1.0)
        @test all(cell.params .<= 1.0)
        @test length(cell.params) == m.cfg["n_cell_params"]
        for conn in cell.conns
            @test conn.input >= -1.0
            @test conn.input <= 1.0
            @test all(conn.state .>= -1.0)
            @test all(conn.state .<= 1.0)
            @test length(conn.state) == m.cfg["n_conn_state"]
            @test all(conn.params .>= -1.0)
            @test all(conn.params .<= 1.0)
            @test length(conn.params) == m.cfg["n_conn_params"]
            @test conn.source[] == cell
            if typeof(conn.dest[]) == Cell
                @test conn.dest[] in m.cells
            else
                in_conns = false
                for c2 in m.cells
                    if conn.dest[] in c2.conns
                        in_conns = true
                        break
                    end
                end
                @test in_conns
            end
        end
    end
    @test nconns <= m.cfg["conns_max"]
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
    nhidden = rand(5:10)
    nreward = rand(5:10)

    @testset "Initialize" begin
        random_init!(m, nin, nout; nhidden=nhidden, nreward=nreward)
        test_model(m)
        @test length(m.cells) == nin + nout + nhidden + nreward
        @test length(m.inputs) == nin
        @test length(m.outputs) == nout
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
        @test all(outputs .>= -1.0)
        @test all(outputs .<= 1.0)
        @test length(outputs) == nout
    end

    @testset "Reward" begin
        reward!(m, rand(nreward))
        test_model(m)
    end

    @testset "Full step" begin
        for i in 1:5
            println("Base controller step $i")
            @time outputs = step!(m, rand(nin))
            reward!(m, rand(nreward))
            test_model(m)
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
    random_init!(m, nin, nout; nreward=nout)

    @testset "Full step" begin
        for i in 1:5
            println(string(cont_constructor, " step ", i))
            @time outputs = step!(m, rand(nin))
            reward!(m, rand(nout))
            test_model(m)
            @test all(outputs .>= -1.0)
            @test all(outputs .<= 1.0)
            @test length(outputs) == nout
        end
    end
end

@testset "Default controller" begin
    controller_test(x->Controller(x))
end

@testset "Const controller" begin
    controller_test(const_controller)
end

@testset "Random controller" begin
    controller_test(rand_controller)
end

# @testset "Static controller" begin
#     controller_test(static_controller;
#                     cfg = Config(Config(), "cfg/static.yaml"))
# end

# @testset "SNN controller" begin
#     controller_test(snn_controller;
#                     cfg = Config(Config(), "cfg/snn.yaml"))
# end

# @testset "CGP controller" begin
#     controller_test(cgp_controller)
# end
