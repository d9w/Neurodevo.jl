# model types

type Cell
  pos::Vector{Float64}
  params::Vector{Int64}
  ctype::Int64
  ntconc::Float64
  ntin::Float64
end

type Synapse
  c1::Cell # reference
  c2::Cell
  weight::Float64
end

type Model
  morphogens::Array{Float64}
  cells::Array{Cell}
  synapses::Array{Synapse}
  itp::InterpGrid
end

type Controller
  division::Function
  child_type::Function
  child_params::Function
  child_position::Function
  apoptosis::Function
  morphogen_diff::Function
  cell_movement::Function
  synapse_formation::Function
  synapse_weight::Function
  nt_output::Function
  nt_update::Function
end

type CellInputs
  ctype::Int64
  params::Vector{Int64}
  ntconc::Float64
end


