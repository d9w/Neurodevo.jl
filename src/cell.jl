export Cell

mutable struct Cell <: Obj
    input::Float64
    output::Float64
    state::Array{Float64}
    params::Array{Float64}
    conns::Array{Conn}
    interface::Bool # is input or output cell
end

function Cell(cfg::Dict, params::Array{Float64};
              interface::Bool=false)
    Cell(0.0, 0.0, zeros(cfg["n_cell_state"]), params,
         Array{Conn}(undef, 0), interface)
end

function get_params(c::Cell, n::Int64)
    c.params
end
