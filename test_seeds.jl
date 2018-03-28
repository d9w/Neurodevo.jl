using CGP
using Base.Test

include("graph_utils.jl")
include("experiments/neural_functions.jl")
include("experiments/seeds.jl")

vstart = -0.65
vthresh = 0.3
thresh = 0.001

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/functions.yaml")

@testset "LIF" begin
    graph = lif_graph(vstart, vthresh)
    chromo = to_chromo(graph)
    for i in 1:10
        inps = 2*rand(5) - 1.0
        lif_outs = min.(1.0, max.(-1.0, lif(inps, vstart, vthresh)))
        chromo_outs = process(chromo, inps)
        @test abs(lif_outs[1] - chromo_outs[1]) < thresh
        @test abs(lif_outs[2] - chromo_outs[2]) < thresh
        @test abs(lif_outs[5] - chromo_outs[5]) < thresh
    end
end
