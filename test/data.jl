include("../evo/data.jl")

function test_classification(cont_constructor::Function;
                             cfg::Dict=Config())
    c = cont_constructor(cfg)
    m = Model(cfg, c)
    X, Y = get_iris(0)
    nin = size(X, 1)
    nout = length(unique(Y))
    layered_init!(m, nin, nout; nhidden=nin, nreward=1)
    fitness = classify(m, X, Y)
    println(fitness)
    @test fitness[1] >= -1.0
end

# @testset "Random controller" begin
#     test_classification(rand_controller)
# end

@testset "Static controller" begin
    test_classification(static_controller; cfg = Config(Config()))
end

# @testset "SNN controller" begin
#     test_classification(snn_controller;
#                         cfg = Config(Config(), "cfg/snn.yaml"))
# end

# @testset "CGP controller" begin
#     include("../evo/cgp.jl")
#     test_classification(cgp_controller)
# end
