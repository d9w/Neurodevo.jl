using JSON
include("plot.jl")
include("evo/data.jl")

function get_individual(filename::String)
    ind_dict = JSON.parse(open(filename, "r"))
    genes = Array{Array{Float64}}(ind_dict["genes"])
    fitness = Array{Float64}(ind_dict["fitness"])
    NeurodevoInd(genes, fitness)
end

function classify_performance(m::Model, X::Array{Float64}, Y::Array{Int64},
                              dataset::String, seed::Int64)
    total_time = 0
    total_memory = 0
    nout = length(unique(Y))
    df = DataFrame(seed=Int[], dataset=String[], epoch=Int[], index=Int[], label=Int[],
                   guess=Int[], certainty=Float64[], correct=Int[])
    for epoch in 1:m.cfg["n_epochs"]
        labels = zeros(Int64, length(Y))
        correct = 0
        for i in eachindex(Y)
            outputs = zeros(nout)
            for s in 1:m.cfg["n_step"]
                outputs, t, bytes, gctime, memallocs = @timed Neurodevo.step!(m, X[:, i])
            end
            labels[i] = argmax(outputs)
            if labels[i] == Y[i]
                correct += 1
                Neurodevo.reward!(m, [1.0])
            end
            push!(df, [seed, dataset, epoch, i, Y[i], labels[i],
                       outputs[labels[i]], correct])
        end
    end
    df
end

function get_performance(ind::NeurodevoInd)
    cfg = cgp_cfg(Config("cfg/evo.yaml"), ind.genes[1])
    c = cgp_controller(cfg, ind.genes[2:end])
    m = Model(cfg, c)

    X, Y = ind.func(ind.seed)
    nin = size(X, 1)
    nout = length(unique(Y))
    nsample = length(Y)
    layered_init!(m, nin, nout; nreward=1, nhidden=cfg["n_hidden"])
    df = classify_performance(m, X, Y, String(split(repr(ind.func), "_")[2]),
                              ind.seed)
    df[:step] = df[:index] + (df[:epoch] .- 1) * nsample
    df[:accuracy] = df[:correct] ./ df[:index]
    df
end

function plot_performance(res::DataFrame; title="Training accuracy",
                          filename="training.pdf")
    plt = plot(res, x=:index, y=:accuracy, color=:epoch,
               Geom.smooth,
               Guide.xlabel("Step"),
               Guide.ylabel("Accuracy"),
               Scale.color_discrete,
               Guide.title(title))
    draw(PDF(filename, 8inch, 6inch), plt)
end

function plot_all_data(ind::NeurodevoInd, id::String)
    ind.func = get_iris
    res = get_performance(ind)
    plot_performance(res; title="Iris training",
                     filename=string(id, "_iris_training.pdf"))
    # ind.func = get_diabetes
    # res = get_performance(ind)
    # plot_performance(res; title="Diabetes training",
    #                  filename=string(id, "_diabetes_training.pdf"))
    # ind.func = get_glass
    # res = get_performance(ind)
    # plot_performance(res; title="Glass training",
    #                  filename=string(id, "_glass_training.pdf"))
end

function plot_all()
    ind = get_individual("gens/1/1000/0100.dna")
    plot_all_data(ind, "1_1000")
    # ind = get_individual("gens/1/2000/0100.dna")
    # plot_all_data(ind, "1_2000")
    # ind = get_individual("gens/89462/1000/0100.dna")
    # plot_all_data(ind, "89462")
end
