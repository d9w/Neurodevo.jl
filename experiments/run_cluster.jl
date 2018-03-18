using RDatasets
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

Logging.configure(filename="stdp.log", level=INFO)

stdp_labels = stdp_cluster(X, Y, 3; n_hidden=50)

acc = randindex(stdp_labels, Y)
