using ArgParse
using CGP
using YAML
include("experiments/clustering.jl")
include("experiments/datasets.jl")
include("experiments/seeds.jl")

settings = ArgParseSettings()
@add_arg_table settings begin
    "--problem"
    arg_type = String
    default = "iris"
    "--experts"
    arg_type = String
    default = "logs/lif_experts.csv"
    "--logfile"
    arg_type = String
    default = "iris_test.log"
    "--seed"
    arg_type = Int64
    default = 0
end

args = parse_args(settings)
Logging.configure(filename=args["logfile"], level=INFO)

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/functions.yaml")
scfg = YAML.load_file("cfg/stdp.yaml")

X, Y = get_data(args["problem"])

function cluster_acc(c::Chromosome, X::Array{Float64}, Y::Array{Int64},
                     pname::String, fname::String, seed=args["seed"])
    n_cluster = length(unique(Y))
    nfunc = i->process(c, i)
    stdp_labels = stdp_cluster(
        X, Y, n_cluster, nfunc; seed=seed, logfile=args["logfile"], problem=pname,
        fname=fname, train_epochs=scfg["train_epochs"],
        weight_mean=scfg["weight_mean"], weight_std=scfg["weight_std"],
        t_train=scfg["t_train"], t_blank=scfg["t_blank"], fr=scfg["fr"],
        pre_target=scfg["pre_target"], stdp_lr=scfg["stdp_lr"],
        stdp_mu=scfg["stdp_mu"], inhib_weight=scfg["inhib_weight"])
    acc = randindex(stdp_labels, Y)
    acc[1]
end

accs = Array{Float64}(10)
srand(args["seed"])
if args["experts"] == "LIF"
    for i in 1:10
        c = to_chromo(lif_graph(-0.65, 0.3))
        accs[i] = cluster_acc(c, X, Y, pname, string(i), i)
    end
else
    experts = readdlm(args["experts"], ',')
    for i in 1:10
        c = CGP.PCGPChromo(experts[i, :], 5, 5)
        accs[i] = cluster_acc(c, X, Y, "iris", string(i))
    end
end

println(accs)
