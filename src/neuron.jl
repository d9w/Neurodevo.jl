type Axon
  position::Vector{Float64}
  age::Int64
end

type Neuron
  position::Vector{Float64}
  axons::Vector{Axon}
  age::Int64
end

function Neuron(dims::Vector{Int64})
  position = rand(length(dims)).*dims
  axons = [Axon(position, 0)]
  Neuron(position, axons, 0)
end

function action!(neuron::Neuron, morphogens::Array, gradients::Array, dims::Vector{Int64})
  nm = length(morphogens)
  str_actions = [[["t_$i"; "a_$i"] for i=1:nm]...;"quiscience";"division";"apoptosis"]
  new_axons = Vector{Axon}()
  actions = ceil(length(str_actions) .* rand(length(neuron.axons)))
  for i in eachindex(neuron.axons)
    axon = neuron.axons[i]
    action = actions[i]
    axon.age += 1
    # 2 is the number of actions per morphogen, as above
    if action <= 2*nm
      axon.position += (2*(action % 2)-1).*gradients[cld(action,2)]
      # stay in bounds
      axon.position = max(min(axon.position, dims),ones(length(dims)))
    elseif action == 2*nm+2 & length(neuron.axons)+length(new_axons) < constants.axon_max
      push!(new_axons, Axon(axon.position, 0))
    end
  end
  axons = neuron.axons[actions .!= 2*nm+3] #apoptosis
  neuron.axons = [axons[:]; new_axons[:]]
  neuron.age += 1
end
