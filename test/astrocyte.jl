include("../src/astrocyte.jl")
include("../src/inits.jl")

function test_creation()
  m = astrocyte_model()
  c = astrocyte_controller()
  print_with_color(:green, "...[passed]\n")
end

function test_step()
  m = astrocyte_model()
  c = astrocyte_controller()
  step!(m, c)
  print_with_color(:green, "...[passed]\n")
end

