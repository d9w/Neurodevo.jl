export static_controller

# A fixed architecture controller with
# + ReLU neurons
# + no neurogenesis or synaptogenesis
# + Reward-modulated Hebbian learning

function static_gen_cell_state_update(cfg::Dict)
    function static_cell_state_update(inputs::Array{Float64})
        # inputs: cell_params, cell_state, inputs
        cell_state = inputs[cfg["n_cell_params"] .+ (1:cfg["n_cell_state"])]
        einputs = @view inputs[(cfg["n_cell_params"]+cfg["n_cell_state"]+1):end]
        outputs = zeros(cfg["n_channels"])
        outputs[1] = max(einputs[1], 0.0)
        outputs[2] = einputs[2]
        cell_state[1] = 0.9 * cell_state[1] + 0.1 * cell_state[1]
        cell_state[2] = 0.9 * cell_state[2] + 0.1 * einputs[2]
        vcat(cell_state, outputs)
    end
end

function static_gen_cell_param_update(cfg::Dict)
    function static_cell_param_update(inputs::Array{Float64})
        cell_params = inputs[1:cfg["n_cell_params"]]
        cell_state = @view inputs[(cfg["n_cell_params"]+1):end]
        cell_params[4] = cell_state[1]
        cell_params[5] = cell_state[2]
        cell_params
    end
end

function static_gen_conn_state_update(cfg::Dict)
    function static_conn_state_update(inputs::Array{Float64})
        # inputs: conn params, state, inputs
        params = @view inputs[1:cfg["n_conn_params"]]
        state = inputs[cfg["n_conn_params"] .+ (1:cfg["n_conn_state"])]
        einputs = @view inputs[(cfg["n_conn_params"]+cfg["n_conn_state"]+1):end]
        state = 0.9 .* state + 0.1 .* einputs
        outputs = params .* einputs
        vcat(state, outputs)
    end
    static_conn_state_update
end

function static_gen_conn_param_update(cfg::Dict)
    function static_conn_param_update(inputs::Array{Float64})
        # inputs: c1 params, c2 params, conn params, state
        params = inputs[2*cfg["n_cell_params"] .+ (1:cfg["n_conn_params"])]
        state = @view inputs[((2*cfg["n_cell_params"] + cfg["n_conn_params"])
                              .+ (1:cfg["n_conn_state"]))]
        dw = state[2] * cfg["lr"] * (
            (state[1] - cfg["weight_mean"]) * (1.0 - params[1]))
        params[1] = max(-1.0, min(1.0, params[1] + dw))
        params
    end
    static_conn_param_update
end

function static_controller(cfg::Dict)
    cell_division = (x...) -> false
    new_cell_params = (x...) -> zeros(cfg["n_cell_params"])
    cell_death = (x...) -> false
    cell_state_update = static_gen_cell_state_update(cfg)
    cell_param_update = static_gen_cell_param_update(cfg)
    connect = (x...) -> false
    new_conn_params = (x...) -> zeros(cfg["n_channels"])
    disconnect = (x...) -> false
    conn_state_update = static_gen_conn_state_update(cfg)
    conn_param_update = static_gen_conn_param_update(cfg)

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
