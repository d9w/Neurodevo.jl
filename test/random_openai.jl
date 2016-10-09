using PyCall
using YAML
using DataFrames
using Gadfly

@pyimport gym
function run_headless(env_name::AbstractString, n_steps::Int64=1000)
  step = 1
  done = false
  r_tot = 0.0

  env = gym.make(env_name)
  # env[:monitor][:start](string(now(),"-random-", env_name))
  o = env[:reset]()
  while !done && step <= n_steps
    action = env[:action_space][:sample]()
    obs, rew, done, info = env[:step](action)
    r_tot += rew
    step += 1
  end
  # env[:monitor][:close]()
  r_tot, step
end

function run_random(group_name::AbstractString)
  config = YAML.load(open("openai.yaml"))
  rtots = []
  steps = []
  ntrial = []
  envs = []
  for env in config["envs"][group_name]
    for i=1:config["n_trials"]
      rtot, step = run_headless(env, config["n_steps"])
      push!(rtots, rtot)
      push!(steps, step)
      push!(ntrial, i)
    end
    append!(envs, [env for i=1:config["n_trials"]])
  end
  results = DataFrame()
  results[:alg] = ["random" for i=1:length(envs)]
  results[:env] = envs
  results[:rtots] = rtots
  results[:steps] = steps
  results
end

function plot_all()
  config = YAML.load(open("openai.yaml"))
  for group_name in keys(config["envs"])
    results = run_random(group_name)
    p = plot(results, x="env", y="rtots", Geom.boxplot)
    img = PNG(string("randomg_",group_name,".png"), 6inch, 4inch)
    draw(img, p)
  end
end

# TODO:
# run 20 trials, plot the progress (two plots) total reward and num steps with x as iteration count
# error bars / area based on 20 trial distribution
# run two - one with a a model that persists between challenges and one that doesn't
# use the same architecture for both
# for random, two trials doesn't matter

# TODO: need input types because locomotion tasks use 2D input
# abstract Input
# type Discrete <: Input
#   min::Int64
#   max::Int64
# end
# type Box <: Input
#   shape::Tuple{Int64}
#   min::Array{Float64}
#   max::Array{Float64}
# end
# input = Input()
# if haskey(env[:action_space],:n)
#   input = min
#   inshape = (1,1)
#   inranges = [0:env[:action_space][:n]]
# elseif haskey(env[:action_space],:shape)
#   n_ins = *(env[:action_space][:shape]...)
# else
#   error("Action space undefined for this environment")
# end
# n_outs = *(env[:observation_space][:shape]...)
