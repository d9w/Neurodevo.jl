using MNIST
using NGE
using Images
using Distances

global const NIN = 7
global const NOUT = 10
global const T_TRAIN = 100
global const T_REST = 40
global const T_STEP = 10
global const B_REWARD = 1.0
const xs = linspace(1.0, DIMS[1], NIN)
const ys = linspace(1.0, DIMS[2], NIN)
const outs = linspace(1.0, DIMS[1], NOUT)

# input the image
function input!(model::Model, cont::Controller, image::Array{Float64})
  # info("input")
  # add value to ntin
  for cell in filter(x->x.params[1] == 1 && x.ctype == 1, model.cells)
    cell.ntin = image[findfirst(xs, cell.pos[1]), findfirst(ys, cell.pos[2])]
  end
end

# return the probability corresponding to each number
function output(model::Model)
  # info("output")
  outputs = zeros(10)
  for cell in filter(x->x.params[1] == 4 && x.ctype == 1, model.cells)
    outputs[findfirst(outs, cell.pos[1])] = cell.ntconc
  end
  outputs /= sum(outputs)
end

# inject morphogen 4 into the final layer of the cell based on correctness
function reward!(model::Model, label::Int64, classification::Array{Float64})
  # info("reward")
  correct = zeros(Float64, 10)
  correct[label] = 1.0
  classification /= sum(classification)
  error = map(x->(correct[x] - classification[x])^2, 1:10)
  for cell in filter(x->x.params[1] == 4 && x.ctype == 1, model.cells)
    ind = findfirst(outs, cell.pos[1])
    reward = 1.0 - error[ind]
    for x in 1:DIMS[1]
      for y in 1:DIMS[2]
        dist = exp(-B_REWARD * euclidean(cell.pos, [x, y, DIMS[3]]))
        model.morphogens[x, y, DIMS[3], 4] += reward * dist
      end
    end
  end
end

# run with astrocytes and morphogen update and supervisory signal
function run!(model::Model, cont::Controller, train::Bool)
  if train
    features, labels = MNIST.traindata(NIN, NIN)
  else
    features, labels = MNIST.testdata(NIN, NIN)
  end
  features /= 256.0
  accuracy = zeros(length(labels))
  for i in eachindex(labels)
    # info("I: $i")
    image = reshape(features[:,i], NIN, NIN)
    label = Integer(labels[i]) + 1
    classification = zeros(Float64, 10)
    # info("training")
    for t = 1:T_TRAIN
      input!(model, cont, image)
      fire!(model, cont)
      classification += output(model)
      if train && (t % T_STEP == 0)
        reward!(model, label, classification)
        step!(model, cont)
      end
    end
    # info("resting")
    for t = 1:T_REST
      fire!(model, cont)
      if train && (t % T_STEP == 0)
        step!(model, cont)
      end
    end
    if indmax(classification) == label
      accuracy[i] = 1.0
    end
    acc = sum(accuracy)
    rate = sum(accuracy) / (i * 1.0)
    info("I: $i A: $acc R: $rate")
  end
  accuracy
end

function train!(model::Model, cont::Controller)
  run!(model, cont, true)
end

function test!(model::Model, cont::Controller)
  run!(model, cont, false)
end

function main()
  info("Creating model and controller")
  model = Model()
  cont = Controller()

  # step once without output for initialization
  info("Initialization")
  step!(model, cont)

  info("Training")
  train_accuracy = train!(model, cont)

  info("Testing")
  test_accuracy = test!(model, cont)
  # TODO: run without supervisory signal
  model, train_accuracy, test_accuracy
end
