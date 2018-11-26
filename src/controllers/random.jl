export rand_controller

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
