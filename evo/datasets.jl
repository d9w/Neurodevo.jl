using DataFrames
using RDatasets

function col_normalize(X::Array{Float64})
    xmin = minimum(X, dims=1)
    xmax = maximum(X, dims=1)
    for i in 1:size(X, 2)
        X[:, i] = (X[:, i] .- xmin[i]) ./ (xmax[i] - xmin[i])
    end
    Array{Float64}(X')
end

function make_x_y(data::DataFrame, seed::Int64)
    x_data = convert(Array{Float64}, data[1:(size(data, 2) - 1)])
    X = col_normalize(x_data)

    labels = unique(data[end])
    Y = convert(Array{Int64}, indexin(data[end], labels))

    rng = MersenneTwister(seed)
    inds = randperm(rng, length(Y))
    X[:, inds], Y[inds]
end

function get_iris(seed::Int64)
    make_x_y(dataset("datasets", "iris"), seed)
end

function get_diabetes(seed::Int64)
    make_x_y(vcat(dataset("MASS", "Pima.te"), dataset("MASS", "Pima.tr")),
             seed)
end

function get_glass(seed::Int64)
    make_x_y(dataset("MASS", "fgl"), seed)
end
