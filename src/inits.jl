export random_init!, layered_init!
# Initializations for models

function random_init!(m::Model, nin::Int64, nout::Int64;
                      nhidden::Int64=nin, nreward::Int64=0)
    nc = m.cfg["n_cell_params"]
    for i in 1:nin
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]); interface=true)
        cell.params[nc] = 0.25
        push!(m.cells, cell)
        push!(m.inputs, Ref(cell))
    end
    for i in 1:nout
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]); interface=true)
        cell.params[nc] = 0.5
        push!(m.cells, cell)
        push!(m.outputs, Ref(cell))
    end
    for i in 1:nreward
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]); interface=true)
        cell.params[nc] = 0.75
        push!(m.cells, cell)
        push!(m.rewards, Ref(cell))
    end
    for i in 1:nhidden
        cell = Cell(m.cfg, rand(m.cfg["n_cell_params"]))
        cell.params[nc] = 1.0
        push!(m.cells, cell)
    end
end

function layer_xy(i::Int64, N::Int64)
    N = ceil(Int64, sqrt(N))
    x = fld(i - 1, N) + 1
    y = i - ((x - 1) * N)
    x / (N + 1), y / (N + 1)
end

function layered_init!(m::Model, nin::Int64, nout::Int64;
                       nhidden::Int64=nin, nreward::Int64=1)
    for i in 1:nin
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nin)
        params[3] = 0.0
        params[end] = 0.25
        cell = Cell(m.cfg, params; interface=true)
        push!(m.cells, cell)
        push!(m.inputs, Ref(cell))
    end
    for i in 1:nhidden
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nhidden)
        params[3] = 0.5
        params[end] = 0.5
        cell = Cell(m.cfg, params)
        push!(m.cells, cell)
    end
    for i in 1:nout
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nout)
        params[3] = 1.0
        params[end] = 0.75
        cell = Cell(m.cfg, params; interface=true)
        push!(m.cells, cell)
        push!(m.outputs, Ref(cell))
    end
    for i in 1:nreward
        params = zeros(m.cfg["n_cell_params"])
        params[1], params[2] = layer_xy(i, nreward)
        params[3] = 0.0
        params[end] = 1.0
        cell = Cell(m.cfg, params; interface=true)
        push!(m.cells, cell)
        push!(m.rewards, Ref(cell))
    end
    all_new_conns = Array{Array{Conn, 1}}(undef, 0)
    for source in m.cells
        new_conns = Array{Conn}(undef, 0)
        for dest in m.cells
            if source.params[3] < dest.params[3] && dest.params[end] != 1.0
                inputs = vcat(source.params, get_params(dest, m.cfg["n_cell_params"]))
                params = m.cont.new_conn_params(inputs)
                conn = Conn(m.cfg, source, dest, params)
                push!(new_conns, conn)
            end
        end
        push!(all_new_conns, new_conns)
    end
    for i in eachindex(m.cells)
        append!(m.cells[i].conns, all_new_conns[i])
    end
    reward = m.rewards[1][]
    for cell in m.cells
        if cell != reward
            for dest in cell.conns
                inputs = vcat(reward.params, get_params(dest, m.cfg["n_cell_params"]))
                params = m.cont.new_conn_params(inputs)
                conn = Conn(m.cfg, reward, dest, params)
                push!(reward.conns, conn)
            end
        end
    end
end
