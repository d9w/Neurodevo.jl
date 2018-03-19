using Clustering
using Logging

struct STDPConfig
    n_hidden::Int64
    dt::Float64
    t_train::Int64
    t_blank::Int64
    fr::Float64
    pre_dt::Float64
    pre_inc::Float64
    pre_target::Float64
    vstart::Float64
    vthresh::Float64
    vmin::Float64
    vmax::Float64
    wmax::Float64
    stdp_lr::Float64
    stdp_mu::Float64
    inhib_weight::Float64
    neuron_function::Function
end

include("network.jl")

function lif(nstate::Array{Float64}, vstart::Float64, vthresh::Float64)
    # Logging.info(@sprintf("N: %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f",
    #                       nstate[1], nstate[2], nstate[3], nstate[4],
    #                       nstate[5], nstate[6]))
    outputs = copy(nstate[1:4])
    if nstate[3] == 0.0
        outputs[1] = vstart
        outputs[2] = 1.0
        outputs[3] = 1.0
    else
        if nstate[6] == 1.0
            outputs[1] = vstart
        else
            outputs[1] = 0.95 * nstate[1] + nstate[5]
        end

        if (outputs[1] >= vthresh) && nstate[2] > 0.8
            outputs[1] = vstart
            outputs[2] = 1.0
        else
            outputs[2] = 0.95 * nstate[2]
        end
    end

    outputs
end


"""
STDP as a clustering function. X is unlabeled data and n is the number of
clusters. X must be between 0 and 1. Returns an Array{Int64} of labels with
size(X,2)
"""
function stdp_cluster(X::Array{Float64}, Y::Array{Int64}, n_cluster::Int64;
                      seed=0, logfile="stdp.log", train_epochs=1,
                      n_hidden=size(X, 1), dt=0.001, t_train=350, t_blank=150,
                      fr=65.0, pre_dt=20.0, pre_inc=1.0, pre_target=0.4,
                      vstart=-65.0, vthresh=30.0, vmin=-100.0, vmax=100.0,
                      wmax=1.0, stdp_lr=0.0001, stdp_mu=2.0,
                      inhib_weight=0.1)::Array{Int64}

    vstart = (vstart - vmin) / (vmax - vmin)
    vthresh = (vthresh - vmin) / (vmax - vmin)

    nfunc = i->lif(i, vstart, vthresh)

    cfg = STDPConfig(n_hidden, dt, t_train, t_blank, fr, pre_dt, pre_inc,
                     pre_target, vstart, vthresh, vmin, vmax, wmax, stdp_lr,
                     stdp_mu, inhib_weight, nfunc)

    srand(seed)
    n_input = size(X, 1)
    network = Network(n_input, n_hidden, n_cluster, cfg)

    # training
    for epoch in 1:train_epochs
        labels = iterate!(network, X, cfg, true)
        acc = randindex(Y, labels)
        Logging.info(@sprintf("R: %d %d %0.4f %0.4f %0.4f %0.4f",
                              epoch, seed, acc[1], acc[2], acc[3], acc[4]))
    end

    # testing
    final_labels = iterate!(network, X, cfg, false)
    acc = randindex(Y, final_labels)
    Logging.info(@sprintf("T: %d %d %0.4f %0.4f %0.4f %0.4f",
                          0, seed, acc[1], acc[2], acc[3], acc[4]))

    final_labels
end
