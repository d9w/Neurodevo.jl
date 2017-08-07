using Base.Test
using E4L

@testset "Individual tests" begin

    N = rand(1:10)
    i = Individual(N, N, N, N)

    new_params = 0.001*randn(N)
    new_state = 0.001*randn(N)

    @testset "Individual creation" begin
        @test i.cell_params == zeros(N, 1)
        @test i.cell_state == zeros(N, 1)
        @test length(i.synapse_params) == 0
        @test length(i.synapse_state) == 0
        @test length(i.synapse_index) == 0
        @test length(i.fitness_history) == 0
        @test length(i.fitness_values) == 0
    end

    @testset "Cell division" begin
        i.controller.cell_division = (x...)->true
        i.controller.child_parameters = (x...)->new_params
        i.controller.child_state = (x...)->new_state
        for i=1:2
            step!(i)
            @test size(i.cell_params) == (N, 2)
            @test size(i.cell_state) == (N, 2)
            @test i.cell_params[:,2] == new_params
            @test i.cell_state[:,2] == new_state
            i.controller.cell_division = (x...)->false
        end
    end

    @testset "Synapse formation" begin
        i.controller.synapse_formation = (x...)->true
        i.controller.synapse_parameters = (x...)->new_params
        i.controller.synapse_state = (x...)->new_state
        for i=1:2
            step!(i)
            @test size(i.synapse_params) == (N, 4)
            @test size(i.synapse_state) == (N, 4)
            @test size(i.synapse_index) == (2, 4)
            @test i.synapse_params == repmat(new_params, 1, 4)
            @test i.synapse_state == repmat(new_state, 1, 4)
            @test i.synapse_index[:,1] == [1,1]
            @test i.synapse_index[:,2] == [1,2]
            @test i.synapse_index[:,3] == [2,1]
            @test i.synapse_index[:,4] == [2,2]
        end
    end

    @testset "Cell update" begin
        i.controller.cell_update = (x...)->zeros(N)
        prev_state = deepcopy(i.cell_state)
        step!(i)
        @test i.cell_state == prev_state + repmat(
            new_state, 1, size(i.cell_state, 2))
    end

    @testset "Synapse update" begin
        i.controller.synapse_update = (x...)->new_state
        prev_state = deepcopy(i.synapse_state)
        step!(i)
        @test i.synapse_state == prev_state + repmat(
            new_state, 1, size(i.cell_state, 2))
    end

    @testset "Input" begin
        i.controller.global_input = (x...)->(1, -1, -1)
        prev_state = deepcopy(i.cell_state)
        inp = rand(20)
        input!(i, inp)
        @test i.cell_state[1, 1] == prev_state[1, 1] + inp[1]
        @test i.cell_state[1, 2] == prev_state[1, 2] + inp[1]
    end

    @testset "Output" begin
        i.controller.global_output = (x...)->(1, -1, 1)
        outs = output(i, 10)
        @test length(outs) == 10
        @test outs[1] == i.cell_state[N, 1] + i.cell_state[N, 2]
        @test outs[2:end] .== 0.0
    end

    @testset "Reward" begin
        i.controller.global_reward = (x...)->(1, 1, -1)
        prev_state = deepcopy(i.cell_state)
        rew = rand(5)
        input!(i, rew)
        @test i.cell_state[1, 1] == prev_state[1, 1] + rew[5]
        @test i.cell_state[1, 2] == prev_state[1, 2] + rew[5]
    end

    @testset "Synapse death" begin
        i.controller.synapse_formation = (x...)->false
        i.controller.synapse_death = (x...)->true
        step!(i)
        @test size(i.sypanse_state, 2) == 0
        @test size(i.sypanse_params, 2) == 0
        @test size(i.sypanse_index, 2) == 0
    end

    @testset "Cell death" begin
        i.controller.synapse_formation = (x...)->true
        i.controller.synapse_death = (x...)->false
        step!(i)
        @test size(i.sypanse_state, 2) == 4
        i.controller.cell_death = (x...)->true
        step!(i)
        @test size(i.cell_state, 2) == 0
        @test size(i.cell_params, 2) == 0
        @test size(i.sypanse_state, 2) == 0
        @test size(i.sypanse_params, 2) == 0
        @test size(i.sypanse_index, 2) == 0
    end
end
