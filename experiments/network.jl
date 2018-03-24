struct Layer
    neurons::Array{Float64}
    weights::Array{Float64}
    inhib::Array{Float64}
    trace::Array{Float64}
end

function Layer(nin::Int64, nout::Int64, cfg::STDPConfig, rng::MersenneTwister)
    neurons = zeros(4, nout)
    weights = min.(1.0, max.(0.0, cfg.weight_std.*randn(rng, nin, nout)
                             .+ cfg.weight_mean))
    inhib = cfg.inhib_weight * (1.0 - eye(nout, nout))
    trace = zeros(nin, nout)
    Layer(neurons, weights, inhib, trace)
end

struct Network
    n_input::Int64
    n_hidden::Int64
    n_output::Int64
    cfg::STDPConfig
    layers::Array{Layer}
end

function Network(n_input::Int64, n_hidden::Int64, n_output::Int64,
                 cfg::STDPConfig; rng::MersenneTwister=MersenneTwister(0))
    hidden = Layer(n_input, n_hidden, cfg, rng)
    output = Layer(n_hidden, n_output, cfg, rng)
    Network(n_input, n_hidden, n_output, cfg, [hidden, output])
end

function spike!(network::Network, layer::Layer, input::BitArray)
    spikes = layer.neurons[1, :] .>= network.cfg.vthresh
    inputs = (input' * layer.weights) - (spikes' * layer.inhib)
    inputs = max.(-1.0, min.(1.0, inputs))
    float_spikes = Float64.(spikes) * 2.0 - 1.0
    for n in 1:size(layer.neurons, 2)
        nstate = [layer.neurons[:, n]; inputs[n]; float_spikes[n]]
        layer.neurons[:, n] = network.cfg.neuron_function(nstate)
    end
    layer.neurons[layer.neurons .< -1.0] .= -1.0
    layer.neurons[layer.neurons .> -1.0] .= 1.0
    spikes
end

function learn!(network::Network, layer::Layer, input_spikes::BitArray,
                output_spikes::BitArray)
    cfg = network.cfg
    layer.trace[input_spikes, :] .+= cfg.pre_inc
    layer.trace .-= (layer.trace ./ cfg.pre_dt)
    layer.weights[:, output_spikes] .+= cfg.stdp_lr .* (
        layer.trace[:, output_spikes] .- cfg.pre_target) .* (
            (cfg.wmax .- layer.weights[:, output_spikes]).^cfg.stdp_mu)
    layer.weights[layer.weights .< 0] .= 0.0
    nothing
end

function step!(network::Network, input_spikes::BitArray, train::Bool)
    output_spikes = copy(input_spikes)
    for l in eachindex(network.layers)
        output_spikes = spike!(network, network.layers[l], input_spikes)
        if train
            learn!(network, network.layers[l], input_spikes, output_spikes)
        end
        input_spikes = copy(output_spikes)
    end
    output_spikes
end

function iterate!(network::Network, X::Array{Float64}, cfg::STDPConfig,
                  train::Bool; rng::MersenneTwister=MersenneTwister(0))
    labels = zeros(Int64, size(X, 2))
    for x in eachindex(labels)
        xfr = X[:, x] * cfg.fr
        out_spikes = zeros(network.n_output)
        input_spikes = rand(rng, network.n_input, cfg.t_train) .< (cfg.dt * xfr)
        for t in 1:cfg.t_train
            out_spikes += step!(network, input_spikes[:, t], train)
        end
        labels[x] = indmax(out_spikes)
        # blank input
        input_spikes = falses(network.n_input)
        for t in 1:cfg.t_blank
            step!(network, input_spikes, train)
        end
    end
    labels
end
