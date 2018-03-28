using CGP

include("graph_utils.jl")

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")

@testset "Chromosome match" begin
    nin = rand(2:5)
    nout = rand(2:5)
    c1 = CGP.PCGPChromo(nin, nout)
    g1 = to_graph(c1)
    c2 = to_chromo(g1)
    @test c1.nin == c2.nin
    @test c1.nout == c2.nout
    @test length(c1.genes) == length(c2.genes)
    @test length(c1.nodes) == length(c2.nodes)
    @test length(c1.outputs) == length(c2.outputs)
    for i in eachindex(c1.outputs)
        @test c1.nodes[c1.outputs[i]].f == c2.nodes[c2.outputs[i]].f
        @test c1.nodes[c1.outputs[i]].active == c2.nodes[c2.outputs[i]].active
        @test c1.nodes[c1.outputs[i]].p == c2.nodes[c2.outputs[i]].p
    end
    c1actives = find([n.active for n in c1.nodes])
    c2actives = find([n.active for n in c2.nodes])
    @test length(c1actives) == length(c2actives)
    for i in 1:10
        ins = rand(c1.nin)
        outs1 = process(c1, ins)
        outs2 = process(c2, ins)
        @test outs1 == outs2
    end
end
@testset "Function match" begin
    func = ins->[*(ins[3], (+(ins[1], ins[2])/2.0)), abs(-(ins[1], ins[3])/2.0),
                 +(ins[1], ins[2])/2.0]
    nin = 3; nout = 3
    fg = get_graph(nin, nout, 3)
    add_node!(fg, 4, CGP.Config.f_sum, x=1, y=2)
    add_node!(fg, 5, CGP.Config.f_mult, x=3, y=4)
    add_node!(fg, 6, CGP.Config.f_aminus, x=1, y=3)
    set_outputs!(fg, nin, nout, 3, [5, 6, 4])
    chromo = to_chromo(fg)
    fg2 = to_graph(chromo)
    @test nv(fg) == nv(fg2)
    @test ne(fg) == ne(fg2)
    for i in 1:10
        ins = rand(chromo.nin)
        outs1 = func(ins)
        outs2 = process(chromo, ins)
        @test outs1 == outs2
    end
end
