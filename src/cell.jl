export Cell

mutable struct Cell
    inputs::Array{Float64}
    outputs::Array{Float64}
    state::Array{Float64}
    params::Array{Float64}
    conns::Array{Conn}
    interface::Bool # is input or output cell
end

function Cell(cfg::Dict, params::Array{Float64};
              interface::Bool=false)
    Cell(zeros(cfg["n_channels"]), zeros(cfg["n_channels"]),
         zeros(cfg["n_cell_state"]), params,
         Array{Conn}(undef, 0), interface)
end
