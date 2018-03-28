function lif(nstate::Array{Float64}, vstart::Float64, vthresh::Float64)
    outputs = zeros(5)
    if (nstate[1] >= vthresh) && nstate[2] < 0.8
        outputs[1] = vstart - nstate[1]
        outputs[2] = 1.0 - nstate[2]
        outputs[5] = 1.0
    else
        outputs[1] = -0.05 * nstate[1] + nstate[5]
        outputs[2] = -0.05 * nstate[2]
    end

    outputs
end

# function izhikevich(nstate::Array{Float64}, vstart::Float64, vthresh::Float64,
#                     vscale::Float64; a::Float64=0.02, b::Float64=0.2,
#                     c::Float64=-65.0, d::Float64=2.0)
#     v = nstate[1] * vscale
#     u = nstate[3] * vscale
#     fired = 0.0
#     if nstate[1] >= vthresh
#         v = c
#         u = u + d
#         fired = 1.0
#     end

#     140 / vscale
#     0.04 * vscale

#     I = nstate[5] * vscale
#     dv = 0.04 * v^2 + 5 * v + 140 - u + I
#     v = v + 0.5 * (0.04 * v^2 + 5 * v + 140 - u + I)
#     v = v + 0.5 * (0.04 * v^2 + 5 * v + 140 - u + I)
#     u = u + a * (b * v - u)

#     [v/vscale-nstate[1], 0.0, u/vscale - nstate[3], 0.0, fired]
# end
