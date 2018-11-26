export Model, set_input!, get_output, reward!

struct Model
    cfg::Dict
    cont::Controller
    cells::Array{Cell}
    inputs::Array{Ref{Cell}}
    outputs::Array{Ref{Cell}}
    state::Array{Int64}
end

function Model(cfg::Dict, cont::Controller)
    cells = Array{Cell}(undef, 0)
    inputs = Array{Ref{Cell}}(undef, 0)
    outputs = Array{Ref{Cell}}(undef, 0)
    state = zeros(1)
    Model(cfg, cont, cells, inputs, outputs, state)
end

function set_input!(m::Model, inputs::Array{Float64})
    for i in eachindex(inputs)
        m.inputs[i][].inputs[1] = inputs[i]
    end
    nothing
end

function get_output(m::Model)
    map(i->m.outputs[i][].outputs[1], eachindex(m.outputs))
end

function reward!(m::Model, reward::Array{Float64})
    for i in eachindex(reward)
        m.outputs[i][].inputs[2] = reward[i]
    end
    nothing
end

function get_conns(m::Model, cell::Cell)
    conns = Array{Conn{Cell}}(undef, 0)
    for c in m.cells
        for conn in c.conns
            if conn.dest[] == cell
                push!(conns, conn)
            end
        end
    end
    conns
end
