include("cgp.jl")
include("darwin.jl")
include("data.jl")

function evaluation(ind::NeurodevoInd)
    X, Y = get_iris()
    nin = size(X, 1)
    nout = length(unique(Y))

    cfg = Config(["cfg/evo.yaml", "cfg/snn.yaml"])
    c = cgp_controller(cfg, ind.genes; cinds=[1 2 3 6 7 8])
    m = Model(cfg, c)
    layered_init!(m, nin, nout)
    for i in 1:m.cfg["T_devo"]
        Neurodevo.step!(m)
    end
    classify(m, X, Y, ind.fitness)
end
