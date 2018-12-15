export step!, develop!

function update_cell_states!(m::Model)
    for cell in m.cells
        inputs = vcat(cell.params, cell.state, cell.input)
        outputs = m.cont.cell_state_update(inputs)
        cell.state = outputs[1:m.cfg["n_cell_state"]]
        cell.output = outputs[end]
    end
end

function update_conn_states!(m::Model)
    for cell in m.cells
        for conn in cell.conns
            inputs = vcat(conn.params, conn.state, cell.output, conn.input)
            outputs = m.cont.conn_state_update(inputs)
            conn.state = outputs[1:m.cfg["n_conn_state"]]
            conn.output = outputs[end]
        end
    end
end

function conn_communication!(m::Model)
    for cell in m.cells
        cell.input = 0
        for conn in cell.conns
            conn.input = 0
        end
    end
    for cell in m.cells
        for conn in cell.conns
            conn.dest[].input += conn.output
        end
    end
    for cell in m.cells
        cell.input = min(1.0, max(-1.0, cell.input))
        for conn in cell.conns
            conn.input = min(1.0, max(-1.0, conn.input))
        end
    end
end

function update_params!(m::Model)
    for cell in m.cells
        inputs = vcat(cell.params, cell.state)
        cell.params = m.cont.cell_param_update(inputs)
        for conn in cell.conns
            inputs = vcat(cell.params, get_params(conn.dest[], m.cfg["n_cell_params"]),
                          conn.params, conn.state)
            conn.params = m.cont.conn_param_update(inputs)
        end
    end
end

function add_cells!(m::Model)
    new_cells = Array{Cell}(undef, 0)
    for cell in m.cells
        div = m.cont.cell_division(cell.params)
        if div && ((length(m.cells) + length(new_cells)) < m.cfg["cells_max"])
            params = m.cont.new_cell_params(cell.params)
            cell = Cell(m.cfg, params)
            push!(new_cells, cell)
        end
    end
    append!(m.cells, new_cells)
    nothing
end

function remove_cells!(m::Model)
    deleted = false
    apop_inds = findall(map(c->~c.interface && m.cont.cell_death(c.params),
                            m.cells))
    if length(apop_inds) > 0
        deleted = true
        deleteat!(m.cells, apop_inds)
    end
    deleted
end

function connect!(m::Model, source::Cell, dest::Obj, new_conns::Array{Conn})
    inputs = vcat(source.params, get_params(dest, m.cfg["n_cell_params"]))
    do_connect = m.cont.connect(inputs)
    if do_connect
        params = m.cont.new_conn_params(inputs)
        conn = Conn(m.cfg, source, dest, params)
        push!(new_conns, conn)
    end
    Int64(do_connect)
end

function connect_cells!(m::Model)
    nconns = 0
    for c in m.cells
        nconns += length(c.conns)
    end
    if nconns == m.cfg["conns_max"]; return nothing; end
    all_new_conns = Array{Array{Conn, 1}}(undef, 0)
    for source in m.cells
        new_conns = Array{Conn}(undef, 0)
        dests = [c.dest[] for c in source.conns]
        for dest in setdiff(m.cells, dests)
            if nconns < m.cfg["conns_max"]
                nconns += connect!(m, source, dest, new_conns)
            end
            for conn in setdiff(dest.conns, dests)
                if nconns < m.cfg["conns_max"]
                    nconns += connect!(m, source, conn, new_conns)
                end
            end
        end
        push!(all_new_conns, new_conns)
    end
    for i in eachindex(m.cells)
        append!(m.cells[i].conns, all_new_conns[i])
    end
end

function disconnect_cells!(m::Model)
    deleted = false
    for source in m.cells
        dinds = Array{Int64}(undef, 0)
        for i in eachindex(source.conns)
            conn = source.conns[i]
            inputs = vcat(source.params,
                          get_params(conn.dest[], m.cfg["n_cell_params"]),
                          conn.params, conn.state)
            do_disconnect = m.cont.disconnect(inputs)
            if do_disconnect
                push!(dinds, i)
            end
        end
        if length(dinds) > 0
            deleted = true
            deleteat!(source.conns, dinds)
        end
    end
    deleted
end

function clean_connections!(m::Model)
    while true
        deleted = false
        for cell in m.cells
            dinds = Array{Int64}(undef, 0)
            for ci in eachindex(cell.conns)
                conn = cell.conns[ci]
                if typeof(conn.dest[]) == Cell
                    if !(conn.dest[] in m.cells)
                        push!(dinds, ci)
                    end
                else
                    in_conns = false
                    for c2 in m.cells
                        if conn.dest[] in c2.conns
                            in_conns = true
                            break
                        end
                    end
                    if !in_conns
                        push!(dinds, ci)
                    end
                end
            end
            if length(dinds) > 0
                deleted = true
                deleteat!(cell.conns, dinds)
            end
        end
        if !deleted
            break
        end
    end
end

function step!(m::Model)
    t = m.state[1]
    update_cell_states!(m)
    update_conn_states!(m)
    conn_communication!(m)
    if mod(t, m.cfg["T_learn"]) == 0
        update_params!(m)
    end
    if mod(t, m.cfg["T_devo"]) == 0
        add_cells!(m)
        connect_cells!(m)
        d1 = remove_cells!(m)
        d2 = disconnect_cells!(m)
        if (d1 || d2)
            clean_connections!(m)
        end
    end
    nothing
end

function step!(m::Model, inputs::Array{Float64})
    set_input!(m, inputs)
    step!(m)
    get_output(m)
end

function develop!(m::Model)
    for i in 1:m.cfg["init_devo"]
        step!(m)
    end
end
