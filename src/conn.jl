export Conn

struct Conn{T}
    inputs::Array{Float64}
    outputs::Array{Float64}
    state::Array{Float64}
    params::Array{Float64}
    source::Ref{T}
    dest::Ref{T}
end

function Conn{T}(cfg::Dict, source::T, dest::T, params::Array{Float64}) where T
    Conn(zeros(cfg["n_channels"]), zeros(cfg["n_channels"]),
         zeros(cfg["n_conn_state"]), params,
         Ref(source), Ref(dest))
end
