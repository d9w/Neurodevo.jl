function import_individual()
    # take an individual from one world and put it in the other
    # reduce that individuals fitness history to the joint set of functions
    # between the two worlds, fix task references accordingly
end

function next_task(f_history::Array{Eval})
    return new_task::Int
end

function mate(f_history_a::Array{Eval}, f_history_b::Array{Eval})
    # Both individuals will be replaced by their children if True
    return mate::Bool
end

function compete(f_history_a::Array{Eval}, f_history_b::Array{Eval})
    # The worse individual will be replaced by the better individual if True
    return compete::Bool
end

function dominates(f_history_a::Array{Eval}, f_history_b::Array{Eval})
    # Compares individuals, A > B if True
    return dominates::Bool
end

function new_state(state::Array{Float}, action::Array{Float})
    return new_state::Array{Float}
end

function reward(state::Array{Float}, action::Array{Float})
    return reward::Float
end
