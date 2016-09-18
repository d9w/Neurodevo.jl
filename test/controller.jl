using StatsBase
using Gadfly
using Colors

include("../src/controller.jl")
include("../src/constants.jl")

n_samples = 1000
n_trials = 1000

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

funcname = "division"
outlabels = ["division"]
inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];"bias"]
infuns = [[x->10*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand()]
inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Float64(x[10]))
profile_graph(division, funcname, outlabels, inlabels, infuns, inform)

funcname = "child_branch"
outlabels = ["child_branch"]
inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];"bias"]
infuns = [[x->10*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand()]
inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Float64(x[10]))
profile_graph(child_branch, funcname, outlabels, inlabels, infuns, inform)

funcname = "child_type"
outlabels = ["child_type"]
inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];"branch"]
infuns = [[x->10*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand(Bool)]
inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Bool(x[10]))
profile_graph(child_type, funcname, outlabels, inlabels, infuns, inform)

funcname = "child_params"
outlabels = ["p$i" for i=1:4]
inlabels = [["m$m" for m=1:4];"pc";["p$i" for i=1:4];"cc";"bias"]
infuns = [[x->10*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];x->rand(1:4);x->rand()]
inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Int64(x[10]),Float64(x[11]))
profile_graph(child_params, funcname, outlabels, inlabels, infuns, inform)

funcname = "child_position"
outlabels = ["pos$i" for i=1:3]
inlabels = [["m$m" for m=1:4];"c";["p$i" for i=1:4];["d$i" for i=1:3];"bias"]
infuns = [[x->10*rand() for i=1:4];x->rand(1:4);[x->rand(1:4) for i=1:4];[x->DIMS[i] for i=1:3];x->rand()]
inform = x->(Array{Float64}(x[1:4]),Int64(x[5]),Array{Int64}(x[6:9]),Array{Float64}(x[10:12]),Float64(x[13]))
profile_graph(child_position, funcname, outlabels, inlabels, infuns, inform)

funcname = "morphogen_diff"
outlabels = ["morphogen_diff"]
inlabels = ["nm";"m";"c";["p$i" for i=1:4];"dist"]
infuns = [x->4;[x->rand(1:4) for i=1:6];x->*(rand(),DIMS...)]
inform = x->(Int64(x[1]),Int64(x[2]),Int64(x[3]),Array{Int64}(x[4:7]),Float64(x[8]))
profile_graph(morphogen_diff, funcname, outlabels, inlabels, infuns, inform)

funcname = "cell_movement"
outlabels = ["pos$i" for i=1:3]
inlabels = [["m$i" for i=1:4];["g$i$j" for i=1:4,j=1:3][:];"c";["p$i" for i=1:4];["d$i" for i=1:3]]
infuns = [[x->10*rand() for i=1:16];x->rand(1:4);[x->rand(1:4) for i=1:4];[x->DIMS[i] for i=1:3]]
inform = x->(Array{Float64}(x[1:4]),reshape(Array{Float64}(x[5:16]),4,3),Int64(x[17]),Array{Int64}(x[18:21]),
             Array{Float64}(x[22:24]))
profile_graph(cell_movement, funcname, outlabels, inlabels, infuns, inform)

funcname = "synapse_weight"
outlabels = ["weight"]
inlabels = [["sm$i" for i=1:4];["am$i" for i=1:4];["sp$i" for i=1:4];["ap$i" for i=1:4]]
infuns = [[x->10*rand() for i=1:8];[x->rand(1:4) for i=1:8]]
inform = x->(Array{Float64}(x[1:4]),Array{Float64}(x[5:8]),Array{Int64}(x[9:12]),Array{Int64}(x[13:16]))
profile_graph(synapse_weight, funcname, outlabels, inlabels, infuns, inform)

funcname = "reward"
outlabels = ["weight"]
inlabels = [["sm$i" for i=1:4];["am$i" for i=1:4];["sp$i" for i=1:4];["ap$i" for i=1:4]]
infuns = [[x->10*rand() for i=1:8];[x->rand(1:4) for i=1:8]]
inform = x->(Array{Float64}(x[1:4]),Array{Float64}(x[5:8]),Array{Int64}(x[9:12]),Array{Int64}(x[13:16]))
profile_graph(synapse_weight, funcname, outlabels, inlabels, infuns, inform)

# funcname = "synapse_formation"
# outlabels = ["formation"]
# inlabels = ["dist"]
# infuns = [x->*(rand(),DIMS...)]
# inform = x->(Float64(x[1]))
# profile_graph(synapse_formation, funcname, outlabels, inlabels, infuns, inform)
