function cell_division(cell::Cell)
    if cell.state[1] < 0.1
        if cell.params[1] == 0.0
            if cell.params[2] < 0.0 && cell.params[2] > -0.2
                if cell.state[1] == 0.25
                    return True
                end
            elseif cell.params[2] == 0.0
                return True
            elseif cell.params[2] > 0.0 && cell.params[2] < 0.2
                if cell.state[1] == 0.75
                    return True
                end
            end
        elseif cell.params[1] == -0.1
            if cell.params[2] < 0.0 && cell.params[2] > -0.1
                if cell.state[1] == 0.25
                    return True
                end
            elseif cell.params[2] == 0.0
                if cell.state[1] == 0.25 || cell.state[1] == 0.75
                    return True
                end
            elseif cell.params[2] > 0.0 && cell.params[2] < 0.1
                if cell.state[1] == 0.75
                    return True
                end
            end
        end
    end
    return False
end

function child_parameters(parent::Cell)
    if parent.state[1] == 0.0
        return [0.1, 0.0]
    elseif parent.state[2] == 0.25
        return [0.0, -0.1]
    elseif parent.state[2] == 0.5
        return [-0.1, 0.0]
    elseif parent.state[2] == 0.75
        return [0.0, 0.1]
    end
    return [0.0, 0.0]
end

function child_state(parent::Cell, params::Array{Float})
    return [-parent.state[1]]
end

function synapse_formation(c1::Cell, c2:Cell)
    return c1.params[1] == c2.params[1] - 0.1
end
