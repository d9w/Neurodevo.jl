using StatsBase
using Gadfly
using Colors

include("../src/model.jl")

function apoptosis()
  # todo, expand to axon/somas
  m = Model()
  apoptosis!(m.cells, m.cells[3])
  @assert ~(3 in collect(keys(m.cells)))
  @assert ~(3 in collect(keys(m.soma_axons)))
end


function handwritten_rules()
  model = Model()
  for i = 1:10
    step!(model)
  end
end
