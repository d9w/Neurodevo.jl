export random_init!, layered_init!
# Initializations for models

function random_init!(m::Model, nin::Int64, nout::Int64; nhidden::Int64=nin)
    for i in 1:nin
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]); interface=true)
        push!(m.cells, cell)
        push!(m.inputs, Ref(cell))
    end
    for i in 1:nout
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]); interface=true)
        push!(m.cells, cell)
        push!(m.outputs, Ref(cell))
    end
    for i in 1:nhidden
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]))
        push!(m.cells, cell)
    end
end

function layer_xy(i::Int64, N::Int64)
    N = ceil(Int64, sqrt(N))
    x = fld(i - 1, N) + 1
    y = i - ((x - 1) * N)
    x / (N + 1), y / (N + 1)
end

function layered_init!(m::Model, nin::Int64, nout::Int64; nhidden::Int64=nin)
    for i in 1:nin
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nin)
        params[4] = 1.0
        cell = Cell(m.cfg, params; interface=true)
        push!(m.cells, cell)
        push!(m.inputs, Ref(cell))
    end
    for i in 1:nout
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nout)
        params[3] = 1.0
        params[4] = 1.0
        cell = Cell(m.cfg, params; interface=true)
        push!(m.cells, cell)
        push!(m.outputs, Ref(cell))
    end
    for i in 1:nhidden
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nhidden)
        params[3] = 0.5
        params[4] = 1.0
        cell = Cell(m.cfg, params)
        push!(m.cells, cell)
    end
end
