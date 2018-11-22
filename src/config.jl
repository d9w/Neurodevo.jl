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
