using StatsBase
using Gadfly
using Colors

include("../src/controller.jl")
include("../src/constants.jl")

const global n_samples = 100
const global n_trials = 100

function profile(cont_func, funcname, outlabels, inlabels, infuns, inform, profgraph=true)
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

  if profgraph
    profile_graph(funcname, inps, outs, cors, inlabels, outlabels)
  end
end

function profile_graph(funcname, inps, outs, cors, inlabels, outlabels)
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

function rand_cell()
  labels = ["id";"p_id";"type";["p$i" for i=1:4];"nt"]
  infuns = [x->rand();x->rand();x->rand(1:N_CTYPES);[x->rand(1:4) for i=1:4];x->rand()]
  inform = x->CellInputs(Float64(x[1]), Float64(x[2]), Int64(x[3]), Array{Int64}(x[4:3+N_PARAMS]),
                         Float64(x[4+N_PARAMS]))
  labels, infuns, inform
end

function test_division(cont::Controller)
  funcname = "division"
  outlabels = ["division"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$m" for m=1:N_MORPHS];cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:N_MORPHS];cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),cell_inform(x[N_MORPHS+1:N_MORPHS+N_PARAMS+4]))
  profile(cont.division, funcname, outlabels, inlabels, infuns, inform)
end

function test_child_type(cont::Controller)
  funcname = "child_type"
  outlabels = ["child_type"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$m" for m=1:N_MORPHS];cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:N_MORPHS];cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),cell_inform(x[N_MORPHS+1:N_MORPHS+N_PARAMS+4]))
  profile(cont.child_type, funcname, outlabels, inlabels, infuns, inform)
end

function test_child_params(cont::Controller)
  funcname = "child_params"
  outlabels = ["p$i" for i=1:4]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$m" for m=1:N_MORPHS];"ccell";cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:N_MORPHS];x->rand(1:N_CTYPES);cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),Int64(x[N_MORPHS+1]),cell_inform(x[N_MORPHS+2:N_MORPHS+N_PARAMS+5]))
  profile(cont.child_params, funcname, outlabels, inlabels, infuns, inform)
end

function test_child_position(cont::Controller)
  funcname = "child_position"
  outlabels = ["pos$i" for i=1:3]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$m" for m=1:N_MORPHS];cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:N_MORPHS];cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),cell_inform(x[N_MORPHS+1:N_MORPHS+N_PARAMS+4]))
  profile(cont.child_position, funcname, outlabels, inlabels, infuns, inform)
end

function test_apoptosis(cont::Controller)
  funcname = "apoptosis"
  outlabels = ["apoptosis"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$m" for m=1:N_MORPHS];cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:N_MORPHS];cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),cell_inform(x[N_MORPHS+1:N_MORPHS+N_PARAMS+4]))
  profile(cont.apoptosis, funcname, outlabels, inlabels, infuns, inform)
end

function test_morphogen_diff(cont::Controller)
  funcname = "morphogen_diff"
  outlabels = ["m$i" for i=1:N_MORPHS]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = ["dist";cell_labels[:]]
  infuns = [x->mean(DIMS)*rand();cell_infuns[:]]
  inform = x->(Float64(x[1]),cell_inform(x[2:N_PARAMS+5]))
  profile(cont.morphogen_diff, funcname, outlabels, inlabels, infuns, inform)
end

function test_cell_movement(cont::Controller)
  funcname = "cell_movement"
  outlabels = ["pos$i" for i=1:N_D]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = [["m$i" for i=1:N_MORPHS];["g$i$j" for i=1:N_MORPHS,j=1:N_D][:];cell_labels[:]]
  infuns = [[x->mean(DIMS)*rand() for i=1:(N_MORPHS*(N_D+1))];cell_infuns[:]]
  inform = x->(Array{Float64}(x[1:N_MORPHS]),reshape(Array{Float64}(x[N_MORPHS+1:N_MORPHS*(N_D+1)]),N_MORPHS,N_D),
               cell_inform(x[N_MORPHS*(N_D+1)+1:N_MORPHS*(N_D+1)+N_PARAMS+4]))
  profile(cont.cell_movement, funcname, outlabels, inlabels, infuns, inform)
end

function test_synapse_fomation(cont::Controller)
  funcname = "synapse_formation"
  outlabels = ["formation"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = ["dist";cell_labels[:];cell_labels[:]]
  infuns = [x->mean(DIMS)*rand();cell_infuns[:];cell_infuns[:]]
  inform = x->(Float64(x[1]),cell_inform(x[2:N_PARAMS+5]),cell_inform(x[N_PARAMS+6:2*N_PARAMS+9]))
  profile(cont.synapse_formation, funcname, outlabels, inlabels, infuns, inform)
end

function test_synapse_weight(cont::Controller)
  funcname = "synapse_weight"
  outlabels = ["weight"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = ["reward";["am$i" for i=1:N_MORPHS];["bm$i" for i=1:N_MORPHS];cell_labels[:];cell_labels[:]]
  infuns = [x->rand();[x->mean(DIMS)*rand() for i=1:2*N_MORPHS];cell_infuns[:];cell_infuns[:]]
  inform = x->(Float64(x[1]),Array{Float64}(x[2:N_MORPHS+1]),Array{Float64}(x[N_MORPHS+2:2*N_MORPHS+1]),
               cell_inform(x[2*N_MORPHS+2:2*N_MORPHS+N_PARAMS+5]),
               cell_inform(x[2*N_MORPHS+N_PARAMS+6:2*N_MORPHS+2*N_PARAMS+9]))
  profile(cont.synapse_weight, funcname, outlabels, inlabels, infuns, inform)
end

function test_synapse_output(cont::Controller)
  funcname = "synapse_output"
  outlabels = ["output"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = ["input";"weight";cell_labels[:]]
  infuns = [x->5*rand();x->randn();cell_infuns[:]]
  inform = x->(Float64(x[1]),Float64(x[2]),cell_inform(x[3:N_PARAMS+6]))
  profile(cont.synapse_output, funcname, outlabels, inlabels, infuns, inform)
end

function test_nt_update(cont::Controller)
  funcname = "nt_update"
  outlabels = ["update"]
  cell_labels, cell_infuns, cell_inform = rand_cell()
  inlabels = ["input";"output";cell_labels[:]]
  infuns = [x->5*rand();x->randn();cell_infuns[:]]
  inform = x->(Float64(x[1]),Float64(x[2]),cell_inform(x[3:N_PARAMS+6]))
  profile(cont.nt_update, funcname, outlabels, inlabels, infuns, inform)
end

function plot_all()
  c = Controller()
  test_division(c)
  test_child_type(c)
  test_child_params(c)
  test_child_position(c)
  test_apoptosis(c)
  test_morphogen_diff(c)
  test_cell_movement(c)
  # test_synapse_fomation(c)
  test_synapse_weight(c)
  test_synapse_output(c)
  test_nt_update(c)
end
