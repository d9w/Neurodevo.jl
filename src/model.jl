export Model, step!, set_input!, get_output, reward!

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

function step!(m::Model)
    nothing
end

function set_input!(m::Model, inputs::Array{Float64})
    nothing
end

function get_output(m::Model)
    rand(length(m.outputs))
end

function step!(m::Model, inputs::Array{Float64})
    set_input!(m, inputs)
    step!(m)
    get_output(m)
end

function reward!(m::Model, reward::Array{Float64})
    nothing
end
