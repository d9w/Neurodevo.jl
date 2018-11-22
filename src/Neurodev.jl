module Neurodev

using Distances
using Random

include("config.jl")
include("controller.jl")
include("controllers/static.jl")
include("controllers/snn.jl")

end
