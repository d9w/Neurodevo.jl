using Test
using Neurodev

@testset "Step tests" begin
    cfg = Config()
    c = Controller(cfg)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)
    random_init!(m, nin, nout)

    @testset "Cell division" begin
        params = rand(m.cfg["n_conn_params"])
        ncells = length(m.cells)
        add_cell!(m, params)
        @test length(m.cells) == ncells + 1
        @test m.cells[end].params == params
        @test length(m.cells[end].conns) == 0
        println("Add cell")
        @time add_cell!(m, params)
    end

    @testset "Connect" begin
        for source in m.cells
            dest = m.cells[rand(eachindex(m.cells))]
            params = rand(m.cfg["n_conn_params"])
            nconns = length(source.conns)
            connect!(m, source, dest, params)
            @test length(source.conns) == nconns + 1
            @test source.conns[end].params == params
            @test source.conns[end].dest[] == dest
        end
    end

    @testset "Cell death" begin
        cell = m.cells[rand(eachindex(m.cells))]
        conns = Neurodev.get_conns(m, cell)
        source = conns[1].source[]
        conns = length(source.conns)
        cells = length(m.cells)
        remove_cell!(m, cell)
        @test length(m.cells) == cells - 1
        @test length(source.conns) == conns - 1
        @test ~(cell in m.cells)
        @test ~(conns[1] in source.conns)
        println("Remove cell")
        cell = m.cells[rand(eachindex(m.cells))]
        @time remove_cell!(m, cell)
    end

    @testset "Disconnect" begin
        source = m.cells[findfirst(
            [length(m.cells[i].conns) > 0 for i in eachindex(m.cells)])]
        conn = source.conns[1]
        nconns = length(source.conns)
        dest = conn.dest[]
        conns = Neurodev.get_conns(m, dest)
        @test conn in conns
        disconnect!(m, conn)
        @test length(source.conns) == nconns - 1
        @test ~(conn in source.conns)
        println("Get conns")
        @time new_conns = Neurodev.get_conns(m, dest)
        @test length(new_conns) == length(conns) - 1
        @test ~(conn in new_conns)
        println("Disconnect")
        source = m.cells[findfirst(
            [length(m.cells[i].conns) > 0 for i in eachindex(m.cells)])]
        conn = source.conns[1]
        @time disconnect!(m, conn)
    end
end
