include("../graph_utils.jl")

function cell_state_update_graph(cfg::Dict)
    nin = Neurodevo.input_lengths(cfg)[4]
    nout = Neurodevo.output_lengths(cfg)[4]
    mg = get_graph(nin, nout, 7)
    add_node!(mg, nin+1, CGP.Config.f_const, p=0.0)
    add_node!(mg, nin+2, CGP.Config.f_max, x=nin+1, y=nin)
    add_node!(mg, nin+3, CGP.Config.f_const, p=0.9)
    add_node!(mg, nin+4, CGP.Config.f_const, p=0.1)
    add_node!(mg, nin+5, CGP.Config.f_mult, x=nin+3, y=cfg["n_cell_params"]+1)
    add_node!(mg, nin+6, CGP.Config.f_mult, x=nin+4, y=nin+2)
    add_node!(mg, nin+7, CGP.Config.f_sum, x=nin+5, y=nin+6)
    set_outputs!(mg, nin, nout, 7, vcat(nin+7, cfg["n_cell_params"] .+ 2:cfg["n_cell_state"]))
    mg
end
