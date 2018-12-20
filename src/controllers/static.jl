export static_controller

# A fixed architecture controller with
# + ReLU neurons
# + no neurogenesis or synaptogenesis
# + Reward-modulated Hebbian learning (Oja's rule)

function static_gen_cell_state_update(cfg::Dict)
    function static_cell_state_update(inputs::Array{Float64})
        # inputs: cell_params, cell_state, input
        cell_state = inputs[cfg["n_cell_params"] .+ (1:cfg["n_cell_state"])]
        einput = inputs[end]
        output = max(einput, 0.0)
        cell_state[1] = 0.9 * cell_state[1] + 0.1 * einput
        vcat(cell_state, output)
    end
end

function static_gen_cell_param_update(cfg::Dict)
    function static_cell_param_update(inputs::Array{Float64})
        # cell params, state -> params
        # save recent cell action in cell params
        cell_params = inputs[1:cfg["n_cell_params"]]
        cell_state = @view inputs[(cfg["n_cell_params"]+1):end]
        cell_params[1] = cell_state[1]
        cell_params
    end
end

function static_gen_new_conn_params(cfg::Dict)
    function new_conn_params(inputs::Array{Float64})
        c1_params = @view inputs[1:cfg["n_conn_params"]]
        c2_params = @view inputs[cfg["n_conn_params"] .+ (1:cfg["n_conn_params"])]
        if c1_params[cfg["n_conn_params"]] == 1.0
            ones(cfg["n_conn_params"])
        else
            rand(cfg["n_conn_params"]) .- 0.5
        end
    end
end

function static_gen_conn_state_update(cfg::Dict)
    function static_conn_state_update(inputs::Array{Float64})
        # inputs: conn params, state, cell in, conn in
        params = @view inputs[1:cfg["n_conn_params"]]
        state = inputs[cfg["n_conn_params"] .+ (1:cfg["n_conn_state"])]
        cell_in = inputs[length(inputs)-1]
        conn_in = inputs[end]
        state[1] = 0.9 * state[1] + 0.1 * conn_in
        output = params[1] * cell_in
        vcat(state, output)
    end
    static_conn_state_update
end

function static_gen_conn_param_update(cfg::Dict)
    function static_conn_param_update(inputs::Array{Float64})
        # inputs: c1 params, c2 params, conn params, state
        c1_params = @view inputs[1:cfg["n_conn_params"]]
        c2_params = @view inputs[cfg["n_conn_params"] .+ (1:cfg["n_conn_params"])]
        params = inputs[2*cfg["n_cell_params"] .+ (1:cfg["n_conn_params"])]
        state = @view inputs[((2*cfg["n_cell_params"] + cfg["n_conn_params"])
                              .+ (1:cfg["n_conn_state"]))]
        dw = (0.01 + state[1]) * c2_params[1] * (c1_params[1] - c2_params[1] * params[1])
        params[1] = max(-1.0, min(1.0, params[1] + dw))
        params
    end
    static_conn_param_update
end

function static_controller(cfg::Dict)
    nouts = output_lengths(cfg)
    cell_division(x::Array{Float64}) = false
    new_cell_params(x::Array{Float64}) = zeros(nouts[2])
    cell_death(x::Array{Float64}) = false
    cell_state_update = static_gen_cell_state_update(cfg)
    cell_param_update = static_gen_cell_param_update(cfg)
    connect(x::Array{Float64}) = false
    new_conn_params = static_gen_new_conn_params(cfg)
    disconnect(x::Array{Float64}) = false
    conn_state_update = static_gen_conn_state_update(cfg)
    conn_param_update = static_gen_conn_param_update(cfg)

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
