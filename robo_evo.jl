using ArgParse
using YAML
using Distributed

@everywhere include("evo/robo.jl")

s = ArgParse.ArgParseSettings()

ArgParse.@add_arg_table(
    s,
    "--seed", arg_type=Int, default=0,
    "--id", arg_type=String, default="test",
    "--log", arg_type=String, default="evolution.log",
    "--cfg", arg_type=String, default="cfg/darwin.yaml",
)

args = ArgParse.parse_args(s)
cfg = YAML.load_file(args["cfg"])
cfg["seed"] = args["seed"]

e = Evolution(NeurodevoInd, cfg; id=args["id"], logfile=args["log"])
e.mutation = uniform_mutation
e.evaluation = robo_eval
e.generation = generation
Darwin.run!(e)
