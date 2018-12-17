include("cgp.jl")
include("darwin.jl")
include("data.jl")

function evaluation(ind::NeurodevoInd)
    cfg = cgp_cfg(Config("cfg/evo.yaml"), ind.genes[1])
    c = cgp_controller(cfg, ind.genes[2:end])
    m = Model(cfg, c)

    X, Y = get_iris()
    nin = size(X, 1)
    nout = length(unique(Y))
    layered_init!(m, nin, nout; nhidden=nin, nreward=1)
    if cfg["init_method"] == 0
        random_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    else
        layered_init!(m, nin, nout; nreward=1, nhidden=cfg["nhidden"])
    end
    develop!(m)
    classify(m, X, Y, ind.fitness)
end
