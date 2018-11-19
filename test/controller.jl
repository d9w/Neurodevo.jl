using Test
using Neurodev

function test_domain(out::Float64)
    @test (out >= -1) && (out <= 1)
end

function test_domain(out::Array{Float64})
    @test all((out .>= -1) .& (out .<= 1))
end

@testset "Controller" begin

    # default configuration
    cfg = Config()

    cell_state = 2.0*rand(cfg.n_cell_state) .- 1.0
    cell_params = 2.0*rand(cfg.n_cell_params) .- 1.0
    branch_state = 2.0*rand(cfg.n_branch_state) .- 1.0
    branch_params = 2.0*rand(cfg.n_branch_params) .- 1.0

    # use default controller for now, switch to test others later
    c = Controller(cfg)

    @testset "Function types" begin
        @test typeof(c.cell_division(
            vcat(cell_params, cell_state))) == Bool
        @test typeof(c.new_cell_params(
            vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.cell_state_update(
            vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.cell_param_update(
            vcat(cell_params, cell_state, rand()))) == Array{Float64,1}
        @test typeof(c.cell_death(
            vcat(cell_params, cell_state))) == Bool

        @test typeof(c.channel_branching(
            vcat(branch_params, branch_state))) == Bool
        @test typeof(c.new_branch_params(
            vcat(branch_params, branch_state))) == Array{Float64,1}
        @test typeof(c.branch_state_update(
            vcat(branch_params, branch_state))) == Array{Float64,1}
        @test typeof(c.branch_param_update(
            vcat(branch_params, branch_state, rand()))) == Array{Float64,1}
        @test typeof(c.branch_connect(
            vcat(branch_params, branch_state, cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.branch_disconnect(
            vcat(branch_params, branch_state, cell_params, cell_state, rand()))) == Bool
        @test typeof(c.branch_pruning(
            vcat(branch_params, branch_state))) == Bool

        @test typeof(c.broadcast_weight(
            vcat(cell_params, cell_params))) == Float64
        @test typeof(c.input(
            vcat(cell_params))) == Array{Float64, 1}
        @test typeof(c.output(
            vcat(cell_params))) == Array{Float64, 1}
    end

    @testset "Outputs sizes" begin
        @test length(c.new_cell_params(
            vcat(cell_params, cell_state))) == cfg.n_cell_params
        @test length(c.cell_state_update(
            vcat(cell_params, cell_state))) == cfg.n_cell_state
        @test length(c.cell_param_update(
            vcat(cell_params, cell_state, rand()))) == cfg.n_cell_params

        @test length(c.new_branch_params(
            vcat(branch_params, branch_state))) == cfg.n_branch_params
        @test length(c.branch_state_update(
            vcat(branch_params, branch_state))) == cfg.n_branch_state
        @test length(c.branch_param_update(
            vcat(branch_params, branch_state, rand()))) == cfg.n_branch_params
        @test length(c.branch_connect(
            vcat(branch_params, branch_state, cell_params, cell_state))) == 2

        @test length(c.input(
            vcat(cell_params))) == 2
        @test length(c.output(
            vcat(cell_params))) == 2
    end

    @testset "Output domain" begin
        test_domain(c.new_cell_params(vcat(cell_params, cell_state)))
        test_domain(c.cell_state_update(vcat(cell_params, cell_state)))
        test_domain(c.cell_param_update(vcat(cell_params, cell_state, rand())))

        test_domain(c.new_branch_params(vcat(branch_params, branch_state)))
        test_domain(c.branch_state_update(vcat(branch_params, branch_state)))
        test_domain(c.branch_param_update(vcat(branch_params, branch_state,
                                               rand())))
        test_domain(c.branch_connect(vcat(branch_params, branch_state,
                                          cell_params, cell_state)))

        test_domain(c.broadcast_weight(vcat(cell_params, cell_params)))
        test_domain(c.input(vcat(cell_params)))
        test_domain(c.output(vcat(cell_params)))
    end
end

nothing
