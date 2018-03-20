function lif(nstate::Array{Float64}, vstart::Float64, vthresh::Float64)
    outputs = nstate[1:4]
    if nstate[3] == 0.0
        outputs[1] = vstart
        outputs[2] = 1.0
        outputs[3] = 1.0
    else
        if nstate[6] == 1.0
            outputs[1] = vstart
        else
            outputs[1] = 0.95 * nstate[1] + nstate[5]
        end

        if (outputs[1] >= vthresh) && nstate[2] > 0.8
            outputs[1] = vstart
            outputs[2] = 1.0
        else
            outputs[2] = 0.95 * nstate[2]
        end
    end

    outputs
end

function izhikevich(nstate::Array{Float64}, vstart::Float64, vthresh::Float64,
                    vscale::Float64; a::Float64=0.02, b::Float64=0.2,
                    c::Float64=-65.0, d::Float64=2.0)
    if nstate[3] == 0.0
        return [c, 0.2 * c, vscale, 0.0] ./ vscale
    end

    v = nstate[1] * vscale
    u = nstate[2] * vscale
    if nstate[6] == 1.0
        v = c
        u = u + d
    end

    I = nstate[5] * vscale
    v = v + 0.5 * (0.04 * v^2 + 5 * v + 140 - u + I)
    v = v + 0.5 * (0.04 * v^2 + 5 * v + 140 - u + I)
    u = u + a * (b * v - u)

    [v, u, vscale, 0.0] ./ vscale
end

function fhn(nstate::Array{Float64}, vstart::Float64, vthresh::Float64;
             scale::Float64=2.0, center::Float64=-0.175,
             a::Float64=0.7, b::Float64=0.8, t::Float64=12.5)
    if nstate[3] == 0.0
        return [vstart, vstart, 1.0, 0.0]
    end

    v = (nstate[1] - center) * 2.0
    u = (nstate[2] - center) * 2.0
    I = (nstate[5] - center) * scale
    v += v - (v.^3) / 3.0 - u + I
    u += 1.0/t * (v + a - b * u)

    [v + center, u + center, scale, 0.0] ./ scale
end
