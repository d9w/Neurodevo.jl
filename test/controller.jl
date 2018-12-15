function test_domain(out::Float64)
    @test (out >= -1) && (out <= 1)
end

function test_domain(out::Array{Float64})
    @test all((out .>= -1) .& (out .<= 1))
end

function scale(ins::Float64)
    2.0 * ins - 1.0
end

function scale(ins::Array{Float64})
    (2.0 .* ins) .- 1.0
end

function test_controller(cfg::Dict, c::Controller)
    cell_state() = scale(rand(cfg["n_cell_state"]))
    cell_params() = scale(rand(cfg["n_cell_params"]))
    conn_state() = scale(rand(cfg["n_conn_state"]))
    conn_params() = scale(rand(cfg["n_conn_params"]))

    @testset "Cell division" begin
        outs = c.cell_division(cell_params())
        @test typeof(outs) == Bool
        @test length(outs) == 1
    end

    @testset "New cell params" begin
        outs = c.new_cell_params(cell_params())
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_cell_params"]
        test_domain(outs)
    end

    @testset "Cell death" begin
        outs = c.cell_death(cell_params())
        @test typeof(outs) == Bool
        @test length(outs) == 1
    end

    @testset "Cell state update" begin
        outs = c.cell_state_update(vcat(cell_params(), cell_state(), rand()))
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_cell_state"] + 1
        test_domain(outs)
    end

    @testset "Cell param update" begin
        outs = c.cell_param_update(vcat(cell_params(), cell_state()))
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_cell_params"]
        test_domain(outs)
    end

    @testset "Connect" begin
        outs = c.connect(vcat(cell_params(), cell_params(), rand()))
        @test typeof(outs) == Bool
        @test length(outs) == 1
    end

    @testset "New conn params" begin
        outs = c.new_conn_params(vcat(cell_params(), cell_params(), rand()))
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_conn_params"]
        test_domain(outs)
    end

    @testset "Disconnect" begin
        outs = c.disconnect(
            vcat(cell_params(), cell_params(), conn_params(), conn_state()))
        @test typeof(outs) == Bool
        @test length(outs) == 1
    end

    @testset "Conn state update" begin
        outs = c.conn_state_update(vcat(conn_params(), conn_state(), rand(2)))
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_conn_state"] + 1
        test_domain(outs)
    end

    @testset "Conn param update" begin
        outs = c.conn_param_update(
            vcat(cell_params(), cell_params(), conn_params(), conn_state()))
        @test typeof(outs) == Array{Float64,1}
        @test length(outs) == cfg["n_conn_params"]
        test_domain(outs)
    end
end

@testset "Default controller" begin
    cfg = Config()
    c = Controller(cfg)
    test_controller(cfg, c)
end

@testset "Random controller" begin
    cfg = Config()
    c = rand_controller(cfg)
    test_controller(cfg, c)
end

@testset "Const controller" begin
    cfg = Config()
    c = const_controller(cfg)
    test_controller(cfg, c)
end

# @testset "Static controller" begin
#     cfg = Config(Config(), "cfg/static.yaml")
#     c = static_controller(cfg)
#     test_controller(cfg, c)
# end

# @testset "SNN controller" begin
#     cfg = Config(Config(), "cfg/snn.yaml")
#     c = snn_controller(cfg)
#     test_controller(cfg, c)
# end

@testset "CGP controller" begin
    include("../evo/cgp.jl")
    cfg = Config()
    c = cgp_controller(cfg)
    test_controller(cfg, c)
end

nothing
