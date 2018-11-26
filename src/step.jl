export add_cell!, remove_cell!, connect!, disconnect!, step!

function add_cell!(m::Model, params::Array{Float64})
    cell = Cell(m.cfg, params)
    push!(m.cells, cell)
    nothing
end

function remove_cell!(m::Model, apop_cell::Cell)
    # delete all connections to the apop cell
    for cell in m.cells
        ind = 0
        for ci in eachindex(cell.conns)
            if cell.conns[ci].dest[] == apop_cell
                ind = ci
                break
            end
        end
        if ind > 0
            deleteat!(cell.conns, ind)
        end
    end
    ind = findfirst(m.cells .== Ref(apop_cell))
    deleteat!(m.cells, ind)
    nothing
end

function connect!(m::Model, source::Cell, dest::Cell, params::Array{Float64})
    conn = Conn{Cell}(m.cfg, source, dest, params)
    push!(source.conns, conn)
    nothing
end

function disconnect!(m::Model, conn::Conn)
    source = conn.source[]
    ind = findfirst(source.conns .== Ref(conn))
    deleteat!(source.conns, ind)
    nothing
end

function step!(m::Model)
    nothing
end

function step!(m::Model, inputs::Array{Float64})
    set_input!(m, inputs)
    step!(m)
    get_output(m)
end
