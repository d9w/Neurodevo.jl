export const_controller

function const_controller(cfg::Dict; c::Float64=1.0)
    cell_division = (x...) -> Bool(round(c))
    new_cell_params = (x...) -> ones(cfg["n_cell_params"]) .* c
    cell_death = (x...) -> Bool(round(c))
    cell_state_update = (x...) -> ones(cfg["n_cell_state"] + cfg["n_channels"]) .* c
    cell_param_update = (x...) -> ones(cfg["n_cell_params"]) .* c
    connect = (x...) -> Bool(round(c))
    new_conn_params = (x...) -> ones(cfg["n_channels"]) .* c
    disconnect = (x...) -> Bool(round(c))
    conn_state_update = (x...) -> ones(cfg["n_conn_state"] + cfg["n_channels"]) .* c
    conn_param_update = (x...) -> ones(cfg["n_conn_params"]) .* c

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
