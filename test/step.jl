@testset "Step tests" begin
    cfg = Config()
    ccon = rand() * 0.49 + 0.51
    c = const_controller(cfg; c=ccon)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)
    random_init!(m, nin, nout)

    @testset "Update cell state" begin
        Neurodevo.update_cell_states!(m)
        for cell in m.cells
            @test all(cell.state .== ccon)
            @test cell.output == ccon
            @test cell.input == 0.0
        end
    end

    @testset "Update conn state" begin
        Neurodevo.update_conn_states!(m)
        for cell in m.cells
            @test cell.input <= 1.0
            @test cell.input >= -1.0
            for conn in cell.conns
                @test all(conn.state .== ccon)
            end
        end
    end

    @testset "Update params" begin
        Neurodevo.update_params!(m)
        for cell in m.cells
            @test all(cell.params .== ccon)
            for conn in cell.conns
                @test all(conn.params .== ccon)
            end
        end
    end

    @testset "Add cells" begin
        ncells = length(m.cells)
        Neurodevo.add_cells!(m)
        @test length(m.cells) == 2 * ncells
        for i in (ncells+1):length(m.cells)
            @test all(m.cells[i].params .== ccon)
        end
    end

    @testset "Remove cells" begin
        Neurodevo.remove_cells!(m)
        @test length(m.cells) == nin + nout
        for i in eachindex(m.cells)
            @test m.cells[i].interface
        end
    end

    @testset "Connect cells" begin
        Neurodevo.connect_cells!(m)
        nconns = 0
        for cell in m.cells
            nconns += length(cell.conns)
        end
        @test nconns <= m.cfg["conns_max"]
    end

    @testset "Disconnect cells" begin
        Neurodevo.disconnect_cells!(m)
        for cell in m.cells
            @test length(cell.conns) == 0
        end
    end
end
