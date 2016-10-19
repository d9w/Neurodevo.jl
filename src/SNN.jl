# SNN.jl
# Spiking Neural Network

global const V_T = -54
global const V_R = -60
global const α = 0.9
global const τp = 13
global const τm = 35
global const τpre = 28
global const τpost = 88
global const Ap = 0.32
global const Am = -0.16

type Neuron
  v::Float64
  v_in::Float64
  spike::Int64
  typ::Int64
end

type Network
  neurons::Vector{Neuron}
  inputs::Vector{Int64}
  outputs::Vector{Int64}
  weights::Array{Float64}
  t::Int64
end

function Network(N_IN::Int64, N_OUT::Int64, N_HIDDEN::Int64)
  N = N_IN + N_OUT + N_HIDDEN
  neurons = Vector{Neuron}(N)
  t = 0
  for i=1:N
    typ = 2
    if i < N_IN
      typ = 0
    elseif i < N_IN+N_OUT
      typ = 1
    end
    neurons[i] = Neuron(0.0, 0.0, 0, t)
  end
  inputs = 1:N_IN
  outputs = N_IN+1:N_OUT
  weights = (V_T-V_R)/2.0*randn(N,N)
  Network(neurons, inputs, outputs, weights, t)
end

function reset(self::Network)
  for n in self.neurons
    n.v = 0
    n.v_in = 0
  end
  nothing
end

function fire(self::Network, ivalues::Vector{Float64})
  @assert length(ivalues) == length(self.inputs)
  for i in eachindex(self.inputs)
    self.neurons[self.inputs[i]].v_in = ivalues[i]
  end

  for n in self.neurons
    n.v = α*(n.v-V_T)*(n.v-V_R) + n.v_in
    n.v_in = 0.0
  end

  for i in eachindex(self.neurons)
    n = self.neurons[i]
    if n.v > V_T
      n.v = V_R
      for j in eachindex(self.neurons)
        self.neurons[j].v_in += self.weights[i,j]
      end
      n.spike = self.t
    end
  end
  self.t += 1
  nothing
end

function weight_update(self::Network, reward::Float64=0.0)
  if reward != 0.0
    ts = Array{Float64}([m.spike-n.spike for n in self.neurons, m in self.neurons])
    eff = Array{Float64}([(1-exp((m.spike-self.t)/τpost))*(1-exp((n.spike-self.t)/τpre))
                          for n in self.neurons, m in self.neurons])
    F = zeros(size(ts))
    F[ts.>0] = Ap*exp(-abs(ts[ts.>0])/τp)
    F[ts.<0] = Am*exp(-abs(ts[ts.<0])/τm)
    self.weights = self.weights .+ (reward * F .* eff)
  end
  nothing
end

function count_output(self::Network, ivalues::Vector{Float64})
  counts = zeros(length(self.outputs))
  for i=1:20
    fire(self, ivalues)
    for o in eachindex(counts)
      if self.neurons[self.outputs[o]].spike == self.t-1
        counts[o] += 1
      end
    end
  end
  findmax(counts)[2]
end

# function plot_reward_response()
#   network = Network(10, 10, 80)
#   firing = [0 0];

#   for t=1:100
#     fire(network, (V_T-V_R)/2.0*rand(10))
#     weight_update(network, (V_T-V_R)/2.0*sin(t/25.0*pi))
#     fired = filter(i->network.neurons[i].spike == network.t-1, eachindex(network.neurons))
#     firing = [firing; [t+0*fired fired]]
#   end

#   infire = firing[(firing[:,2].>0) .* (firing[:,2].<11),:];
#   outfire = firing[(firing[:,2].>10) .* (firing[:,2].<21),:];
#   hiddenfire = firing[(firing[:,2].>20) .* (firing[:,2].<81),:];

#   l = @layout [a; b{0.2h}]
#   p = scatter(infire[:,1], infire[:,2], yticks=0:10:60, markersize=4, markercolor=:red, label="input", layout=l)
#   scatter!(p[1], outfire[:,1], outfire[:,2], yticks=0:10:60, markersize=4, markercolor=:green, label="output")
#   scatter!(p[1], hiddenfire[:,1], hiddenfire[:,2], yticks=0:10:60, markersize=4, markercolor=:black, label="hidden")

#   plot!(p[2], (V_T-V_R)/2.0*sin(1:100/25.0*pi))

#   pyplot(leg=false, ticks=nothing)
#   x = y = linspace(-5, 5, 40)
#   zs = zeros(0,40)

#   @gif for i in linspace(0, 2π, 100)
#       f(x,y) = sin(x + 10sin(i)) + cos(y)

#       # create a plot with 3 subplots and a custom layout
#       l = @layout [a{0.7w} b; c{0.2h}]
#       p = plot(x, y, f, st = [:surface, :contourf], layout=l)

#       # add a tracking line
#       fixed_x = zeros(40)
#       z = map(f,fixed_x,y)
#       plot!(p[1], fixed_x, y, z, line = (:black, 5, 0.2))
#       vline!(p[2], [0], line = (:black, 5))

#       # add to and show the tracked values over time
#       zs = vcat(zs, z')
#       plot!(p[3], zs, alpha = 0.2, palette = cgrad(:blues).colors)
#   end
# end
