export const_controller

function const_controller(cfg::Dict; c::Float64=1.0)
    nouts = output_lengths(cfg)
    cell_division(x::Array{Float64}) = Bool(round(c))
    new_cell_params(x::Array{Float64}) = ones(nouts[2]) .* c
    cell_death(x::Array{Float64}) = Bool(round(c))
    cell_state_update(x::Array{Float64}) = ones(nouts[4]) .* c
    cell_param_update(x::Array{Float64}) = ones(nouts[5]) .* c
    connect(x::Array{Float64}) = Bool(round(c))
    new_conn_params(x::Array{Float64}) = ones(nouts[7]) .* c
    disconnect(x::Array{Float64}) = Bool(round(c))
    conn_state_update(x::Array{Float64}) = ones(nouts[9]) .* c
    conn_param_update(x::Array{Float64}) = ones(nouts[10]) .* c

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
