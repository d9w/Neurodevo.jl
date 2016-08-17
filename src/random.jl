include("model.jl")

g = DiGraph()
#m_decay = 0.0:0.1:0.5
constants = Constants(0.1, 100, 1000, 0.2, 20)

for i=1:30
  N = 10*i
  model = Model([10,10,10], N)
  step_count = 0
  change = 0
  avg_change = 0

  while step_count < constants.step_total && (change > constants.change_stop || step_count < constants.step_init)
    g = model.synapses
    step!(model)
    change = ne(difference(g, model.synapses))/ne(model.synapses)
    avg_change += change
    step_count += 1
  end

  avg_change /= step_count
  stats = evaluate(model)
  print_joined(STDOUT, [N; stats...; step_count; avg_change], ",")
  print("\n")
  save("../graphs/graph$N", g, :lg)
end
