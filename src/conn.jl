export Conn
import Base.show

abstract type Obj end

mutable struct Conn <: Obj
    input::Float64
    output::Float64
    state::Array{Float64}
    params::Array{Float64}
    source::Ref
    dest::Ref
end

function Conn(cfg::Dict, source::Obj, dest::Obj, params::Array{Float64})
    Conn(0.0, 0.0, zeros(cfg["n_conn_state"]), params, Ref(source), Ref(dest))
end

function get_params(c::Conn, n::Int64)
    vcat(c.params, zeros(n - length(c.params)))
end

function show(io::IO, x::Conn)
    print(io, string("Conn(", x.input, ", ", x.output, ", ", x.state, ", ",
                     x.params, ", ", x.source[].params, ", ", x.dest[].params, ")"));
end
