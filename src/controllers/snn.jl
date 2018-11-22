export snn_controller

# A fully-connected snn baseline controller with
# + LIF neurons
# + neurogenesis
# + synaptogenesis (on neurogenesis)
# + RM-STDP

function snn_gen_cell_division(cfg::Dict)
    function snn_cell_division(cell_params::Array{Float64})
        ((cell_params[4] > cfg["divide_age"]) &&
         (cell_params[5] > cfg["divide_act"]))
     end
     snn_cell_division
end

function snn_new_cell_params(cell_params::Array{Float64})
    vcat(copy(cell_params[1:3]) + 0.01 * randn(3), 1.0, 0.0)
end

function snn_gen_cell_death(cfg::Dict)
    function snn_cell_death(cell_params::Array{Float64})
        cell_params[5] < cfg["cell_apop"]
    end
    snn_cell_death
end

function snn_gen_cell_state_update(cfg::Dict)
    function snn_cell_state_update(inputs::Array{Float64})
        # inputs: cell_params, cell_state, inputs
        cell_state = inputs[cfg["n_cell_params"] .+ (1:cfg["n_cell_state"])]
        einputs = @view inputs[(cfg["n_cell_params"]+cfg["n_cell_state"]+1):end]
        spike = cell_state[1] > (cfg["v_thresh"] + cell_state[2])
        outputs = zeros(cfg["n_channels"])
        if spike
            cell_state[1] = cfg["v_reset"]
            cell_state[2] = min(1.0, cell_state[2] + cfg["v_thresh_inc"])
            outputs[1] = cfg["spike_exc"]
            outputs[2] = -cfg["spike_inhib"]
        else
            cell_state[1] = min(1.0, max(-1.0,
                cfg["v_leak"] * cell_state[1] + einputs[1] + einputs[2]))
            cell_state[2] = cfg["v_thresh_dec"] * cell_state[2]
        end
        vcat(cell_state, outputs)
    end
    snn_cell_state_update
end

function snn_gen_cell_param_update(cfg::Dict)
    function snn_cell_param_update(inputs::Array{Float64})
        cell_params = copy(inputs[1:cfg["n_cell_params"]])
        cell_state = @view inputs[(cfg["n_cell_params"]+1):end]
        cell_params[4] *= 0.999
        cell_params[5] = cell_state[2]
        cell_params
    end
    snn_cell_param_update
end

function snn_connect(inputs::Array{Float64})
    true
end

function snn_gen_new_conn_params(cfg::Dict)
    function snn_new_conn_params(inputs::Array{Float64})
        # inputs: c1 params, c2 params
        c1 = @view inputs[1:cfg["n_cell_params"]]
        c2 = @view inputs[(cfg["n_cell_params"]+1):end]
        weights = zeros(cfg["n_conn_params"])
        if c2[3] > c1[3]
            weights[1] = randn()*cfg["weight_std"] + cfg["weight_mean"]
        end
        weights[2] = exp(-euclidean(c1[1:3], c2[1:3]))
        weights
    end
    snn_new_conn_params
end

function snn_disconnect(inputs::Array{Float64})
    false
end

function snn_gen_conn_state_update(cfg::Dict)
    function snn_conn_state_update(inputs::Array{Float64})
        # inputs: conn params, state, inputs
        params = @view inputs[1:cfg["n_conn_params"]]
        state = inputs[cfg["n_conn_params"] .+ (1:cfg["n_conn_state"])]
        einputs = @view inputs[(cfg["n_conn_params"]+cfg["n_conn_state"]+1):end]
        state .*= cfg["trace_dec"]
        nstate = min.(1.0, max.(-1.0, einputs + state))
        outputs = params .* einputs
        vcat(nstate, outputs)
    end
    snn_conn_state_update
end

function snn_gen_conn_param_update(cfg::Dict)
    function snn_conn_param_update(inputs::Array{Float64})
        # inputs: c1 params, c2 params, conn params, state
        params = inputs[2*cfg["n_cell_params"] .+ (1:cfg["n_conn_params"])]
        state = @view inputs[((2*cfg["n_cell_params"] + cfg["n_conn_params"])
                              .+ (1:cfg["n_conn_state"]))]
        dw = state[2] * cfg["stdp_lr"] * (
            (state[1] - cfg["pre_target"]) * (1.0 - params[1]))
        params[1] = max(-1.0, min(1.0, params[1] + dw))
        params
    end
    snn_conn_param_update
end

function snn_controller(cfg::Dict)
    Controller(snn_gen_cell_division(cfg),
               snn_new_cell_params,
               snn_gen_cell_death(cfg),
               snn_gen_cell_state_update(cfg),
               snn_gen_cell_param_update(cfg),
               snn_connect,
               snn_gen_new_conn_params(cfg),
               snn_disconnect,
               snn_gen_conn_state_update(cfg),
               snn_gen_conn_param_update(cfg))
end
