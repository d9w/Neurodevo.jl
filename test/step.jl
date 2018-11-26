using Test
using Neurodev

@testset "Step tests" begin
    cfg = Config()
    ccon = rand() * 0.51 + 0.49
    c = const_controller(cfg; c=ccon)
    m = Model(cfg, c)
    nin = rand(5:10)
    nout = rand(5:10)
    random_init!(m, nin, nout)

    @testset "Update cell state" begin
        Neurodev.update_cell_states!(m)
        for cell in m.cells
            @test all(cell.state .== ccon)
            @test all(cell.outputs .== ccon)
            @test all(cell.inputs .== 0.0)
        end
    end

    @testset "Update conn state" begin
        Neurodev.update_conn_states!(m)
        for cell in m.cells
            @test all(cell.inputs .<= 1.0)
            @test all(cell.inputs .>= -1.0)
            for conn in cell.conns
                @test all(conn.state .== ccon)
            end
        end
    end

    @testset "Update params" begin
        Neurodev.update_params!(m)
        for cell in m.cells
            @test all(cell.params .== ccon)
            for conn in cell.conns
                @test all(conn.params .== ccon)
            end
        end
    end

    @testset "Add cells" begin
        ncells = length(m.cells)
        Neurodev.add_cells!(m)
        @test length(m.cells) == 2 * ncells
        for i in (ncells+1):length(m.cells)
            @test all(m.cells[i].params .== ccon)
        end
    end

    @testset "Remove cells" begin
        Neurodev.remove_cells!(m)
        @test length(m.cells) == nin + nout
        for i in eachindex(m.cells)
            @test m.cells[i].interface
        end
    end

    @testset "Connect cells" begin
        Neurodev.connect_cells!(m)
        for cell in m.cells
            @test length(cell.conns) == length(m.cells)
        end
    end

    @testset "Disconnect cells" begin
        Neurodev.disconnect_cells!(m)
        for cell in m.cells
            @test length(cell.conns) == 0
        end
    end
end
