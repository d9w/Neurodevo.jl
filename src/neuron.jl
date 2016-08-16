type Axon
  position::Vector{Int64}
  age::Int64
end

type Neuron
  position::Vector{Int64}
  axons::Vector{Axon}
  connections::Vector{Int64}
  age::Int64
end

function Neuron(dims::Vector{Int64})
  position = map(x -> rand(1:x), dims)
  axons = [Axon(position, 0)]
  Neuron(position, axons, [], 0)
end

function emission(neuron::Neuron, morphogens::Vector{Float64})
  return rand(length(morphogens))
end

function action!(neuron::Neuron, morphogens::Vector{Float64}, dims::Vector{Int64})
  nm = length(morphogens)
  str_actions = [["m$i" for i=1:nm]...,"quiscience","division","apoptosis"]
  actions = ceil(length(str_actions) .* rand(length(neuron.axons)))
  for axon in neuron.axons[actions .<= nm]
    axon.position += round(randn(length(axon.position)))
    axon.position = max(min(axon.position, dims),ones(length(dims)))
  end
  for axon in neuron.axons
    axon.age += 1
  end
  if length(neuron.axons) < 20
    for axon in neuron.axons[actions .== nm+2]
      push!(neuron.axons, Axon(axon.position, 0))
    end
  end
  # axons = neuron.axons[actions .!= nm+3]
  # for axon in axons
  #   axon.age += 1
  # end
  # neuron.axons = [axons[:]; new_axons[:]]
  neuron.age += 1
end
