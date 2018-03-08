using Base.Test
using E4L

@testset "Controller" begin

    n_cell_state = rand(1:10)
    n_cell_params = rand(1:10)
    n_syn_state = rand(1:10)
    n_syn_params = rand(1:10)
    cell_state = 2.0*rand(n_cell_state)-1.0
    cell_params = 2.0*rand(n_cell_params)-1.0
    syn_state = 2.0*rand(n_syn_state)-1.0
    syn_params = 2.0*rand(n_syn_params)-1.0

    # use default controller for now, switch to test others later
    c = Controller(n_cell_state, n_cell_params, n_syn_state, n_syn_params)

    @testset "Function types" begin
        @test typeof(c.cell_division(vcat(cell_params, cell_state))) == Bool
        @test typeof(c.child_parameters(vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.child_state(vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.cell_update(vcat(
            cell_params, cell_state, syn_params, syn_state))) == Array{Float64,1}
        @test typeof(c.cell_death(vcat(cell_params, cell_state))) == Bool
        @test typeof(c.synapse_formation(vcat(
            cell_params, cell_state, cell_params, cell_state))) == Bool
        @test typeof(c.synapse_parameters(vcat(
            cell_params, cell_state, cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.synapse_state(vcat(
            cell_params, cell_state, cell_params,
            cell_state, syn_params))) == Array{Float64,1}
        @test typeof(c.synapse_update(vcat(
            cell_params, cell_state, cell_params,
            cell_state, syn_params))) == Array{Float64,1}
        @test typeof(c.synapse_death(vcat(
            cell_params, cell_state, cell_params, cell_state, syn_params))) == Bool
        @test typeof(c.input(vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.output(vcat(cell_params, cell_state))) == Array{Float64,1}
        @test typeof(c.reward(vcat(cell_params, cell_state))) == Array{Float64,1}
    end

    @testset "Outputs sizes" begin
        @test length(c.child_parameters(vcat(
            cell_params, cell_state))) == n_cell_params
        @test length(c.child_state(vcat(
            cell_params, cell_state))) == n_cell_state
        @test length(c.cell_update(vcat(
            cell_params, cell_state, syn_params, syn_state))) == n_cell_state
        @test length(c.synapse_parameters(vcat(
            cell_params, cell_state, cell_params, cell_state))) == n_syn_params
        @test length(c.synapse_state(vcat(
            cell_params, cell_state, cell_params, cell_state, syn_params))) == n_syn_state
        @test length(c.synapse_update(vcat(
            cell_params, cell_state, cell_params, cell_state, syn_params))) == n_syn_state
        @test length(c.input(vcat(
            cell_params, cell_state))) == 3
        @test length(c.output(vcat(
            cell_params, cell_state))) == 3
        @test length(c.reward(vcat(
            cell_params, cell_state))) == 3
    end

    @testset "Output domain" begin
        ans = c.child_parameters(vcat(cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.child_state(vcat(cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.cell_update(vcat(cell_params, cell_state, syn_params,syn_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.synapse_parameters(vcat(cell_params, cell_state, cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.synapse_state(vcat(
            cell_params, cell_state, cell_params, cell_state, syn_params))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.synapse_update(vcat(
            cell_params, cell_state, cell_params, cell_state,syn_params))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.input(vcat(cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.output(vcat(cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
        ans = c.reward(vcat(cell_params, cell_state))
        @test all((ans .>= -1) .& (ans .<= 1))
    end
end

nothing
