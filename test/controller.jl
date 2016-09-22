using StatsBase
using Gadfly
using Colors

include("../src/controller.jl")
include("../src/constants.jl")

const global n_samples = 1000
const global n_trials = 1000

function profile_graph(cont_func, funcname, outlabels, inlabels, infuns, inform)
  println(funcname)
  inps = Array{Float64}(n_samples*n_trials*length(inlabels),length(inlabels))
  outs = Array{Float64}(n_samples*n_trials*length(inlabels),length(outlabels))
  cors = Array{Float64}(length(inlabels)*n_trials,length(outlabels))
  point = 1
  pc = 1
  for t=1:n_trials
    ins = [infuns[i](0.0) for i in eachindex(infuns)]
    for i=1:length(inlabels)
      incopy = copy(ins)
      p1 = copy(point)
      for s=1:n_samples
        incopy[i] = infuns[i](0.0)
        out = cont_func(inform(incopy)...)
        if length(out) > 1
          out = vec(out)
        end
        inps[point,:] = incopy
        outs[point,:] = out
        point+=1
      end
      for o=1:length(outlabels)
        cors[pc,o] = abs(corspearman(inps[p1:p1+n_samples-1,i],outs[p1:p1+n_samples-1,o]))
      end
      pc += 1
    end
  end
  cors[isnan(cors)] = 0.0;

  ps = Array{Compose.Context}(length(outlabels))
  for o in eachindex(outlabels)
    pv = plot(x=ones(n_samples), y=outs[:,o], Geom.violin, Guide.xlabel(""), Guide.xticks(label=false))
    pc = plot(x=repmat(inlabels,n_trials), y=cors[:,o], Geom.histogram2d(xbincount=length(inlabels),ybincount=10),
              Scale.color_continuous(colormap=p->RGB(0.0,0.75-0.75*p,1.0-p)),
              Guide.xlabel(outlabels[o]), Guide.xticks(orientation=:vertical),
              Guide.yticks(ticks=collect(0.0:0.2:1.0)), Guide.colorkey(""))
    if o==1
      append!(pv.guides,[Guide.ylabel("Value")])
      append!(pc.guides,[Guide.ylabel("Correlation")])
    else
      append!(pv.guides,[Guide.ylabel("")])
      append!(pc.guides,[Guide.ylabel("")])
    end
    pc.theme.key_position=:right
    ps[o] = vstack(pv, pc)
  end
  draw(PNG("/home/d9w/Documents/projects/axon-guidance/plots/$funcname.png",
          (10*length(outlabels))cm, 20cm), hstack(ps...))
end

function plot_division(cont::Controller)
  funcname = "division"
  outlabels = ["division"]
  inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];"velo"]
  infuns = [[x->mean(DIMS)*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand()]
  inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Float64(x[10]))
  profile_graph(cont.division, funcname, outlabels, inlabels, infuns, inform)
end

function plot_child_type(cont::Controller)
  funcname = "child_type"
  outlabels = ["child_type"]
  inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4]]
  infuns = [[x->mean(DIMS)*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4]]
  inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]))
  profile_graph(cont.child_type, funcname, outlabels, inlabels, infuns, inform)
end

function plot_child_params(cont::Controller)
  funcname = "child_params"
  outlabels = ["p$i" for i=1:4]
  inlabels = [["m$m" for m=1:4];"pc";["p$i" for i=1:4];"cc"]
  infuns = [[x->mean(DIMS)*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand(1:4)]
  inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Int64(x[10]))
  profile_graph(cont.child_params, funcname, outlabels, inlabels, infuns, inform)
end

function plot_child_position(cont::Controller)
  funcname = "child_position"
  outlabels = ["pos$i" for i=1:3]
  inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4]]
  infuns = [[x->mean(DIMS)*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4]]
  inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]))
  profile_graph(cont.child_position, funcname, outlabels, inlabels, infuns, inform)
end

function plot_apoptosis(cont::Controller)
  funcname = "apoptosis"
  outlabels = ["apoptosis"]
  inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];"velo"]
  infuns = [[x->mean(DIMS)*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand()]
  inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Float64(x[10]))
  profile_graph(cont.apoptosis, funcname, outlabels, inlabels, infuns, inform)
end

function plot_morphogen_diff(cont::Controller)
  funcname = "morphogen_diff"
  outlabels = ["morphogen_diff"]
  inlabels = ["nm";"m";"c";["p$i" for i=1:4];"dist"]
  infuns = [x->4;[x->rand(1:4) for i=1:6];x->*(rand(),DIMS...)]
  inform = x->(Int64(x[1]),Int64(x[2]),Int64(x[3]),Array{Int64}(x[4:7]),Float64(x[8]))
  profile_graph(cont.morphogen_diff, funcname, outlabels, inlabels, infuns, inform)
end

function plot_cell_movement(cont::Controller)
  funcname = "cell_movement"
  outlabels = ["pos$i" for i=1:3]
  inlabels = [["m$i" for i=1:4];["g$i$j" for i=1:4,j=1:3][:];"c";["p$i" for i=1:4]]
  infuns = [[x->mean(DIMS)*rand() for i=1:16];x->rand(1:4);[x->rand(1:4) for i=1:4]]
  inform = x->(Array{Float64}(x[1:4]),reshape(Array{Float64}(x[5:16]),4,3),Int64(x[17]),Array{Int64}(x[18:21]))
  profile_graph(cont.cell_movement, funcname, outlabels, inlabels, infuns, inform)
end

function plot_synapse_weight(cont::Controller)
  funcname = "synapse_weight"
  outlabels = ["weight"]
  inlabels = [["sm$i" for i=1:4];["am$i" for i=1:4];["sp$i" for i=1:4];["ap$i" for i=1:4];"reward"]
  infuns = [[x->mean(DIMS)*rand() for i=1:8];[x->rand(1:4) for i=1:8];x->rand()]
  inform = x->(Array{Float64}(x[1:4]),Array{Float64}(x[5:8]),Array{Int64}(x[9:12]),Array{Int64}(x[13:16]),
              Float64(x[17]))
  profile_graph(cont.synapse_weight, funcname, outlabels, inlabels, infuns, inform)
end

function plot_synapse_fomation(cont::Controller)
  funcname = "synapse_formation"
  outlabels = ["formation"]
  inlabels = ["dist"]
  infuns = [x->*(rand(),DIMS...)]
  inform = x->(Float64(x[1]))
  profile_graph(cont.synapse_formation, funcname, outlabels, inlabels, infuns, inform)
end

function plot_all()
  c = Controller()
  plot_division(c)
  plot_child_type(c)
  plot_child_params(c)
  plot_child_position(c)
  plot_apoptosis(c)
  plot_morphogen_diff(c)
  plot_cell_movement(c)
  plot_synapse_weight(c)
  # plot_synapse_fomation(c)
end
