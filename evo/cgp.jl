using Neurodevo
using CGP

CGP.Config.init("cfg/cgp/base.yaml")
CGP.Config.init("cfg/cgp/functions.yaml")

function make_chromos!(cfg::Dict; cinds::Array{Int64}=collect(1:10))
    chromos = Array{CGP.PCGPChromo}(undef, length(cinds))
    nins = Neurodevo.input_lengths(cfg)
    nouts = Neurodevo.output_lengths(cfg)
    for i in eachindex(chromos)
        chromos[i] = CGP.PCGPChromo(nins[cinds[i]], nouts[cinds[i]])
    end
    cfg["chromosomes"] = chromos
    chromos
end

function make_chromos!(cfg::Dict, genes::Array{Array{Float64}};
                       cinds::Array{Int64}=collect(1:10))
    chromos = Array{CGP.PCGPChromo}(undef, length(cinds))
    nins = Neurodevo.input_lengths(cfg)
    nouts = Neurodevo.output_lengths(cfg)
    for i in eachindex(cinds)
        chromos[i] = CGP.PCGPChromo(genes[i], nins[cinds[i]], nouts[cinds[i]])
    end
    cfg["chromosomes"] = chromos
    chromos
end

function make_controller(cfg::Dict)
    chromos = cfg["chromosomes"]
    snnc = snn_controller(cfg)
    cell_division(x::Array{Float64}) = process(chromos[1], x)[1] > 0
    new_cell_params(x::Array{Float64}) = process(chromos[2], x)
    cell_death(x::Array{Float64}) = process(chromos[3], x)[1] > 0
    # cell_state_update(x::Array{Float64}) = process(chromos[1], x)
    # cell_param_update(x::Array{Float64}) = process(chromos[2], x)
    connect(x::Array{Float64}) = process(chromos[4], x)[1] > 0
    # connect(x::Array{Float64}) = true
    new_conn_params(x::Array{Float64}) = process(chromos[5], x)
    # new_conn_params(x::Array{Float64}) = rand(cfg["n_cell_params"])
    disconnect(x::Array{Float64}) = process(chromos[6], x)[1] > 0
    # conn_state_update(x::Array{Float64}) = process(chromos[3], x)
    # conn_param_update(x::Array{Float64}) = process(chromos[4], x)

    Controller(cell_division, new_cell_params, cell_death,
               snnc.cell_state_update, snnc.cell_param_update,
               connect, new_conn_params, disconnect,
               snnc.conn_state_update, snnc.conn_param_update)
    # Controller(cell_division, new_cell_params, cell_death,
    #            cell_state_update, cell_param_update,
    #            connect, new_conn_params, disconnect,
    #            conn_state_update, conn_param_update)
end

function cgp_controller(cfg::Dict; cinds::Array{Int64}=collect(1:10))
    make_chromos!(cfg; cinds=cinds)
    make_controller(cfg)
end

function cgp_controller(cfg::Dict, genes::Array{Array{Float64}};
                        cinds::Array{Int64}=collect(1:10))
    make_chromos!(cfg, genes; cinds=cinds)
    make_controller(cfg)
end
