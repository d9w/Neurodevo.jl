include("model.jl")

constants = Constants(0.1);
g = DiGraph()

for i=1:1000
  model = Model()
  step_count = 0
  change = 0

  while step_count < 1000 && (change > .10 || step_count < 100)
    g = model.synapses
    step!(model)
    change = ne(difference(g, model.synapses))/ne(model.synapses)
    step_count += 1
  end

  stats = evaluate(model)
  print_joined(STDOUT, [stats...; step_count], ",")
  print("\n")
  save("../graphs/graph$i", g, :lg)
end
