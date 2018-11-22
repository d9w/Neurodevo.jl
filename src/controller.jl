export Controller, rand_controller

struct Controller
    cell_division::Function # cell params -> bool
    new_cell_params::Function # cell params -> child params
    cell_death::Function # cell params -> bool
    cell_state_update::Function # cell params, state, inputs -> state, outputs
    cell_param_update::Function # cell params, state -> params

    connect::Function # cell 1 params, cell 2 params -> bool
    new_conn_params::Function # cell 1 params, cell 2 params -> conn params
    disconnect::Function # c1 params, c2 params, conn params, state -> bool
    conn_state_update::Function # conn params, state, inputs -> state, outputs
    conn_param_update::Function # c1 params, c2 params, conn params, state -> params
end

function Controller(cfg::Dict)
    cell_division = (x...) -> false
    new_cell_params = (x...) -> zeros(cfg["n_cell_params"])
    cell_death = (x...) -> false
    cell_state_update = (x...) -> zeros(cfg["n_cell_state"] + cfg["n_channels"])
    cell_param_update = (x...) -> zeros(cfg["n_cell_params"])
    connect = (x...) -> false
    new_conn_params = (x...) -> zeros(cfg["n_channels"])
    disconnect = (x...) -> false
    conn_state_update = (x...) -> zeros(cfg["n_conn_state"] + cfg["n_channels"])
    conn_param_update = (x...) -> zeros(cfg["n_conn_params"])

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end

function rand_controller(cfg::Dict)
    Random.seed!(cfg["seed"])
    cell_division = (x...) -> rand(Bool)
    new_cell_params = (x...) -> rand(cfg["n_cell_params"])
    cell_death = (x...) -> rand(Bool)
    cell_state_update = (x...) -> rand(cfg["n_cell_state"] + cfg["n_channels"])
    cell_param_update = (x...) -> rand(cfg["n_cell_params"])
    connect = (x...) -> rand(Bool)
    new_conn_params = (x...) -> rand(cfg["n_channels"])
    disconnect = (x...) -> rand(Bool)
    conn_state_update = (x...) -> rand(cfg["n_conn_state"] + cfg["n_channels"])
    conn_param_update = (x...) -> rand(cfg["n_conn_params"])

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
