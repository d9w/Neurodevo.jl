using ArgParse
using YAML
using Distributed

const IN_SLURM = "SLURM_JOBID" in keys(ENV)
IN_SLURM && using ClusterManagers

if IN_SLURM
    pids = addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])))
else
    pids = addprocs()
end

@everywhere include("evo/data.jl")

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
cfg["n_fitness"] = 10

e = Evolution(NeurodevoInd, cfg; id=args["id"], logfile=args["log"])
e.mutation = uniform_mutation
e.evaluation = data_evaluation
run!(e)
