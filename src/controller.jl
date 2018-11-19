export Controller

struct Controller
    cell_division::Function # cell params, state -> bool
    new_cell_params::Function # parent cell params, state -> child params
    cell_state_update::Function # cell params, state -> new state
    cell_param_update::Function # cell params, state, nbranches -> new params
    cell_death::Function # cell params, state -> bool
    channel_branching::Function # branch params, state -> bool
    new_branch_params::Function # parent branch params, state -> child params
    branch_state_update::Function # branch params, state -> new state
    branch_param_update::Function # branch params, state, nconns -> new params
    branch_connect::Function # branch params, state, cell params, state -> bool
    branch_disconnect::Function # branch params, state, cell params, state, weight -> bool
    branch_pruning::Function # branch params, state, nconns -> bool
    broadcast_weight::Function # cell 1 params, cell 2 params -> weight
    input::Function # cell params -> bool, index
    output::Function # cell params -> bool, index
end

function Controller(cfg::Config)
    cell_division = (x...) -> false
    new_cell_params = (x...) -> zeros(cfg.n_cell_params)
    cell_state_update = (x...) -> zeros(cfg.n_cell_state)
    cell_param_update = (x...) -> zeros(cfg.n_cell_params)
    cell_death = (x...) -> false
    channel_branching = (x...) -> false
    new_branch_params = (x...) -> zeros(cfg.n_branch_params)
    branch_state_update = (x...) -> zeros(cfg.n_branch_state)
    branch_param_update = (x...) -> zeros(cfg.n_branch_params)
    branch_connect = (x...) -> zeros(2)
    branch_disconnect = (x...) -> false
    branch_pruning = (x...) -> false
    broadcast_weight = (x...) -> 0.0
    input = (x...) -> zeros(2)
    output = (x...) -> zeros(2)

    Controller(cell_division, new_cell_params, cell_state_update,
               cell_param_update, cell_death, channel_branching,
               new_branch_params, branch_state_update, branch_param_update,
               branch_connect, branch_disconnect, branch_pruning,
               broadcast_weight, input, output)
end
