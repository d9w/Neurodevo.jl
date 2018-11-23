export Conn

struct Conn{T}
    inputs::Array{Float64}
    outputs::Array{Float64}
    state::Array{Float64}
    params::Array{Float64}
    source::Ref{T}
    cells::Array{Ref{T}}
end

function Conn{T}(cfg::Dict, source::Ref{T}, params::Array{Float64}) where T
    Conn(zeros(cfg["n_channels"]), zeros(cfg["n_channels"]),
         zeros(cfg["n_conn_state"]), params,
         source, Array{Conn{T}})
end
