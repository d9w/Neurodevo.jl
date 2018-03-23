using ArgParse
using CGP
using YAML
include("clustering.jl")
include("datasets.jl")

settings = ArgParseSettings()
@add_arg_table settings begin
    "--seed"
    arg_type = Int
    default = 0
    "--logfile"
    arg_type = String
    default = "stdp.log"
end

args = parse_args(settings)
srand(args["seed"])
Logging.configure(filename=args["logfile"], level=INFO)

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")
scfg = YAML.load_file("cfg/stdp.yaml")
fname = string(args["seed"])

problems = ["iris", "spirals"]
data = Dict()
for p in problems
    X, Y = get_data(p)
    dp = Dict()
    dp[:X] = X
    dp[:Y] = Y
    data[p] = dp
end

test_problem = "yeast"
test_X, test_Y = get_data("yeast")

function cluster_acc(c::Chromosome, X::Array{Float64}, Y::Array{Int64},
                     pname::String)
    n_cluster = length(unique(Y))
    nfunc = i->process(c, i)
    stdp_labels = stdp_cluster(
        X, Y, n_cluster, nfunc; seed=0, logfile=args["logfile"], problem=pname,
        fname=fname, train_epochs=scfg["train_epochs"],
        weight_mean=scfg["weight_mean"], weight_std=scfg["weight_std"],
        t_train=scfg["t_train"], t_blank=scfg["t_blank"], fr=scfg["fr"],
        pre_target=scfg["pre_target"], stdp_lr=scfg["stdp_lr"],
        stdp_mu=scfg["stdp_mu"], inhib_weight=scfg["inhib_weight"])
    acc = randindex(stdp_labels, Y)
    acc[1]
end

function cluster_fit(c::Chromosome)
    fit = 0.0
    for p in problems
        fit += cluster_acc(c, data[p][:X], data[p][:Y], p)
    end
    fit /= length(problems)
end

function gen_fit(c::Chromosome)
    cluster_acc(c, test_X, test_Y, test_problem)
end

maxfit, best = oneplus(CGPChromo, 6, 4, cluster_fit; seed=args["seed"],
                       record_best=true, record_fitness=gen_fit)
Logging.info(@sprintf("E%0.6f", -maxfit))
