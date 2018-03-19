using RDatasets
using ArgParse
include("clustering.jl")

iris = dataset("datasets", "iris")

iris_data = Array{Float64}([iris[:SepalLength] iris[:SepalWidth] iris[:PetalWidth] iris[:PetalLength]]);
xmin = minimum(iris_data, 1)
xmax = maximum(iris_data, 1)
X = Array{Float64}(size(iris_data, 2), size(iris_data, 1))
for i in 1:size(iris_data, 2)
    X[i, :] = (iris_data[:, i] .- xmin[i]) ./ (xmax[i] - xmin[i])
end

sps = unique(iris[:Species])
iris[:label] = indexin(iris[:Species], sps)
Y = Array{Int64}(iris[:label])

settings = ArgParseSettings()
@add_arg_table settings begin
    "--seed"
    arg_type = Int
    default = 0
    "--logfile"
    arg_type = String
    default = "stdp.log"
    "--train_epochs"
    arg_type = Int
    default = 1
    "--n_hidden"
    arg_type = Int
    default = 40
    "--t_train"
    arg_type = Int
    default = 350
    "--t_blank"
    arg_type = Int
    default = 150
    "--fr"
    arg_type = Float64
    default = 65.0
    "--pre_target"
    arg_type = Float64
    default = 0.4
    "--stdp_lr"
    arg_type = Float64
    default = 0.0001
    "--stdp_mu"
    arg_type = Float64
    default = 2.0
    "--inhib_weight"
    arg_type = Float64
    default = 0.1
end
settings

args = parse_args(settings)
Logging.configure(filename=args["logfile"], level=INFO)

stdp_labels = stdp_cluster(
    X, Y, 3; seed=args["seed"], logfile=args["logfile"],
    train_epochs=args["train_epochs"], n_hidden=args["n_hidden"],
    t_train=args["t_train"], t_blank=args["t_blank"], fr=args["fr"],
    pre_target=args["pre_target"], stdp_lr=args["stdp_lr"],
    stdp_mu=args["stdp_mu"], inhib_weight=args["inhib_weight"])

acc = randindex(stdp_labels, Y)

Logging.info(@sprintf("E%0.6f", -acc[1]))
