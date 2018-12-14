include("../evo/data.jl")

function test_classification(cont_constructor::Function; cfg::Dict=Config())
    c = cont_constructor(cfg)
    m = Model(cfg, c)
    X, Y = get_iris()
    nin = size(X, 1)
    nout = length(unique(Y))
    layered_init!(m, nin, nout; nhidden=1, nreward=0)
    for i in 1:m.cfg["T_devo"]
        step!(m)
    end

    pfits = zeros(5)
    fitness = classify(m, X, Y, pfits)
    println(fitness)
    @test fitness[1] >= 0.0
    @test length(fitness) == length(pfits)
end

@testset "Random controller" begin
    test_classification(rand_controller)
end

# @testset "Static controller" begin
#     test_classification(static_controller;
#                         cfg = Config(Config(), "cfg/static.yaml"))
# end

# @testset "SNN controller" begin
#     test_classification(snn_controller;
#                         cfg = Config(Config(), "cfg/snn.yaml"))
# end

# @testset "CGP controller" begin
#     include("../evo/cgp.jl")
#     test_classification(cgp_controller)
# end
