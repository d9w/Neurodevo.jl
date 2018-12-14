export Config
using YAML

function Config(;filename="cfg/base.yaml")
    YAML.load_file(filename)
end

function Config(filename::String)
    YAML.load_file(filename)
end

function Config(filenames::Array{String})
    merge([YAML.load_file(f) for f in filenames]...)
end

function Config(cfg::Dict)
    cfg
end

function Config(cfg::Dict, other::Dict)
    merge(cfg, other)
end

function Config(cfg::Dict, other::String)
    merge(cfg, Config(other))
end

function Config(cfg::Dict, other::Array{String})
    merge(cfg, Config(other))
end

function input_lengths(cfg::Dict)
    inputs = Array{Int64}(undef, 10)
    inputs[1] = cfg["n_cell_params"]
    inputs[2] = cfg["n_cell_params"]
    inputs[3] = cfg["n_cell_params"]
    inputs[4] = cfg["n_cell_params"] + cfg["n_cell_state"] + 1
    inputs[5] = cfg["n_cell_params"] + cfg["n_cell_state"]
    inputs[6] = 2*cfg["n_cell_params"]
    inputs[7] = 2*cfg["n_cell_params"]
    inputs[8] = 2*cfg["n_cell_params"] + cfg["n_conn_params"] + cfg["n_conn_state"]
    inputs[9] = cfg["n_conn_params"] + cfg["n_conn_state"] + 2
    inputs[10] = 2*cfg["n_cell_params"] + cfg["n_conn_params"] + cfg["n_conn_state"]
    inputs
end

function output_lengths(cfg::Dict)
    outputs = Array{Int64}(undef, 10)
    outputs[1] = 1
    outputs[2] = cfg["n_cell_params"]
    outputs[3] = 1
    outputs[4] = cfg["n_cell_state"] + 1
    outputs[5] = cfg["n_cell_params"]
    outputs[6] = 1
    outputs[7] = cfg["n_conn_params"]
    outputs[8] = 1
    outputs[9] = cfg["n_conn_state"] + 1
    outputs[10] = cfg["n_conn_params"]
    outputs
end
