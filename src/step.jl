export step!

function update_cell_states!(m::Model)
    for cell in m.cells
        inputs = vcat(cell.params, cell.state, cell.inputs)
        outputs = m.cont.cell_state_update(inputs)
        cell.state = outputs[1:m.cfg["n_cell_state"]]
        cell.outputs = outputs[(m.cfg["n_cell_state"]+1):end]
        cell.inputs .= 0
    end
end

function update_conn_states!(m::Model)
    for cell in m.cells
        for conn in cell.conns
            inputs = vcat(conn.params, conn.state, cell.outputs)
            outputs = m.cont.conn_state_update(inputs)
            conn.state = outputs[1:m.cfg["n_conn_state"]]
            conn.dest[].inputs += outputs[(m.cfg["n_conn_state"]+1):end]
        end
    end
    for cell in m.cells
        cell.inputs = min.(1.0, max.(-1.0, cell.inputs ./ m.cfg["cells_max"]))
    end
end

function update_params!(m::Model)
    for cell in m.cells
        inputs = vcat(cell.params, cell.state)
        cell.params = m.cont.cell_param_update(inputs)
        for conn in cell.conns
            inputs = vcat(cell.params, conn.dest[].params, conn.params,
                          conn.state)
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
    apop_inds = findall(map(c->~c.interface && m.cont.cell_death(c.params),
                            m.cells))
    if length(apop_inds) > 0
        apop_cells = m.cells[apop_inds]
        for cell in m.cells
            inds = findall(map(ci->cell.conns[ci].dest[] in apop_cells,
                               eachindex(cell.conns)))
            if length(inds) > 0
                deleteat!(cell.conns, inds)
            end
        end
        deleteat!(m.cells, apop_inds)
    end
    nothing
end

function connect_cells!(m::Model)
    for source in m.cells
        dests = [c.dest[] for c in source.conns]
        for dest in m.cells
            if ~(dest in dests)
                inputs = vcat(source.params, dest.params)
                do_connect = m.cont.connect(inputs)
                if do_connect
                    params = m.cont.new_conn_params(inputs)
                    conn = Conn{Cell}(m.cfg, source, dest, params)
                    push!(source.conns, conn)
                end
            end
        end
    end
end

function disconnect_cells!(m::Model)
    for source in m.cells
        dinds = Array{Int64}(undef, 0)
        for i in eachindex(source.conns)
            conn = source.conns[i]
            inputs = vcat(source.params, conn.dest[].params, conn.params,
                          conn.state)
            do_disconnect = m.cont.disconnect(inputs)
            if do_disconnect
                push!(dinds, i)
            end
        end
        deleteat!(source.conns, dinds)
    end
end

function step!(m::Model)
    t = m.state[1]
    update_cell_states!(m)
    update_conn_states!(m)
    if mod(t, m.cfg["T_learn"]) == 0
        update_params!(m)
    end
    if mod(t, m.cfg["T_devo"]) == 0
        add_cells!(m)
        remove_cells!(m)
        connect_cells!(m)
        disconnect_cells!(m)
    end
    nothing
end

function step!(m::Model, inputs::Array{Float64})
    set_input!(m, inputs)
    step!(m)
    get_output(m)
end
