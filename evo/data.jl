using DataFrames
using RDatasets
using Neurodevo
using Random

function col_normalize(X::Array{Float64})
    xmin = minimum(X, dims=1)
    xmax = maximum(X, dims=1)
    for i in 1:size(X, 2)
        X[:, i] = (X[:, i] .- xmin[i]) ./ (xmax[i] - xmin[i])
    end
    Array{Float64}(X')
end

function get_iris()
    iris = dataset("datasets", "iris")

    iris_data = Array{Float64}(hcat(iris[:SepalLength], iris[:SepalWidth],
                                    iris[:PetalWidth], iris[:PetalLength]));
    X = col_normalize(iris_data)

    sps = unique(iris[:Species])
    iris[:label] = indexin(iris[:Species], sps)
    Y = Array{Int64}(iris[:label])
    rng = MersenneTwister(1234)
    inds = randperm(rng, length(Y))
    X[:, inds], Y[inds]
end

function classify(m::Model, X::Array{Float64}, Y::Array{Int64},
                  pfits::Array{Float64})
    base_fit = pfits[1]
    new_fits = zeros(length(pfits))
    total_time = 0
    total_memory = 0
    for epoch in eachindex(pfits)
        labels = zeros(Int64, length(Y))
        for i in eachindex(Y)
            outputs, t, bytes, gctime, memallocs = @timed Neurodevo.step!(m, X[:, i])
            labels[i] = argmax(outputs)
            if labels[i] == Y[i]
                Neurodevo.reward!(m, [1.0])
            end
            total_time += t
            total_memory += bytes
            if (total_time >= m.cfg["time_max"] ||
                total_memory >= m.cfg["memory_max"])
                break
            end
        end
        fit = sum(labels .== Y) / length(Y)
        new_fits[epoch] = fit
        base_fit = fit
        if (fit < pfits[epoch] || fit < base_fit ||
            total_time >= m.cfg["time_max"] ||
            total_memory >= m.cfg["memory_max"])
            break
        end
    end
    new_fits
end

