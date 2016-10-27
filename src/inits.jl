# Initializations for models

include("model.jl")

function random_model(N::Int64=100)
  m = Model()
  for n=1:N
    add_cell!(m, Cell())
  end
  m
end

function add_layer!(model::Model, lnum::Int64, maxl::Int64, nn::Int64)
  n = floor(Int64, sqrt(nn))
  z = lnum/maxl
  for x in linspace(0.0, 1.0, n)
    for y in linspace(0.0, 1.0, n)
      add_cell!(model, Cell([x, y, z], [lnum], 1, 0.0, 0.0))
      add_cell!(model, Cell([x, y, z], [lnum], 2, 0.0, 0.0))
    end
  end
end

function astrocyte_model()
  m = Model()
  add_layer!(m, 1, 4, 7*7)#28*28)
  add_layer!(m, 2, 4, 5*5)#20*20)
  add_layer!(m, 3, 4, 4*4)#10*10)
  add_layer!(m, 4, 4, 9)
  m
end

function permute_model()
  m = Model()
  param_permutes = Array{Int64}(N_MORPHS^N_PARAMS, N_PARAMS)
  i = 1
  for c = Counter([N_MORPHS for m=1:N_PARAMS])
    param_permutes[i,:] = c
    i += 1
  end
  n_cuts = ceil(Int64,N_MORPHS^(N_PARAMS/N_D))
  cuts = [linspace(1,DIMS[d],n_cuts) for d=1:N_D]
  i = 1
  for c = Counter([n_cuts for m=1:N_D])
    if i <= size(param_permutes)[1]
      pos = convert(Array{Float64},[cuts[d][c[d]] for d=1:N_D])
      params = vec(param_permutes[i,:])
      add_cell!(m.cells, Cell(pos, params, 1, 0.0, 0.0))
      i += 1
    else
      break
    end
  end
  m.cells = cells
end
