using Neurodev
using Random

function rand_controller(cfg::Config)
    Random.seed!(cfg.seed)
    cell_division = (x...) -> rand(Bool)
    new_cell_params = (x...) -> rand(cfg.n_cell_params)
    cell_state_update = (x...) -> rand(cfg.n_cell_state)
    cell_param_update = (x...) -> rand(cfg.n_cell_params)
    cell_death = (x...) -> rand(Bool)
    channel_branching = (x...) -> rand(Bool)
    new_branch_params = (x...) -> rand(cfg.n_branch_params)
    branch_state_update = (x...) -> rand(cfg.n_branch_state)
    branch_param_update = (x...) -> rand(cfg.n_branch_params)
    branch_connect = (x...) -> rand(2)
    branch_disconnect = (x...) -> rand(Bool)
    branch_pruning = (x...) -> rand(Bool)
    broadcast_weight = (x...) -> rand()
    input = (x...) -> rand(2)
    output = (x...) -> rand(2)

    Controller(cell_division, new_cell_params, cell_state_update,
               cell_param_update, cell_death, channel_branching,
               new_branch_params, branch_state_update, branch_param_update,
               branch_connect, branch_disconnect, branch_pruning,
               broadcast_weight, input, output)
end
