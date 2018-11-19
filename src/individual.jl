export Individual

struct Individual
    chromosomes::Array{Array{Float64}}
    cfg::Config
    model::Model
    controller::Controller
end

function Individual(chromosomes::Array{Float64})
    cfg = make_config(chromosomes[1])
end
