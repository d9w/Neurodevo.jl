using GLVisualize, GLWindow, GeometryTypes, FileIO
using GLAbstraction, Colors, Reactive, ProgressMeter

include("../src/astrocyte.jl")
include("../src/inits.jl")

window = glscreen()
timesignal = loop(linspace(0f0, 1f0, 360))

cube = HyperRectangle(Vec3f0(0), Vec3f0(0.04))
n = 20
const wx,wy,wz = widths(cube)

mesh = GLNormalMesh(cube)
# mesh = GLNormalMesh(loadasset("/home/d9w/Documents/projects/graphing/GLVisualize.jl/assets/cat.obj"))
timepi = const_lift(*, timesignal, 2f0*pi)

model = astrocyte_model()
cont = astrocyte_controller()

poses = rand(n, 3)
mind = zeros(3)
maxd = ones(3)

# function position(t, x)
#   # poses[x,:] = min(max(vec(poses[x,:])+0.01*randn(3),mind),maxd)
#   pos = Point3f0(((model.cells[x].pos)./DIMS)...)
#   # pos = Point3f0(x*(sqrt(wx^2+wy^2))*abs(cos(x*y*t*pi/360)), -y*wy*abs(cos(x*y*t*pi/180)), y*wz*abs(cos(t*y*pi/90)))
#   # pos = Point3f0(x*wx, sin(x*t)*wy, sin(t)*wz)
#   # dir = Point3f0(0, wy, wz)
#   # pos = pos + sin(t)*dir - 1.0
#   # Point3f0(model.cells[x].pos...)
# end

position_signal = map(timepi) do t
  # step!(model, cont)
  model.morphogens += 0.1*rand(size(model.morphogens)...)
  # vec(Point3f0[Point3f0(model.cells[x].pos...) for x in keys(model.cells)])
  vec(Point3f0[Point3f0(((cell.pos)./DIMS)...) for cell in model.cells])
  # vec(Point3f0[position(t,x) for x=1:n])
end

# function scale_gen(v0, nv)
# 	l = length(v0)
# 	@inbounds for i=1:l
# 		v0[i] = Vec3f0(1,1,sin((nv*l)/i))/2
# 	end
# 	v0
# end
# function color_gen(v0, t)
# 	l = length(v0)
# 	@inbounds for x=1:l
# 		# v0[x] = RGBA{U8}(x/l,(cos(t)+1)/2,(sin(x/l/3)+1)/2.,1.)
# 		v0[x] = RGBA{U8}(0.35,0.422,(sin(x/l/3)+1)/2.,1.)
# 	end
# 	v0
# end

color_signal = map(timepi) do t
  maxmorph = maximum(model.morphogens)
  vec([RGBA{U8}(([model.itp[c.pos...,m] for m=1:3]./maxmorph)...,0.8) for c in model.cells])
end

t = const_lift(x->x+0.1, timesignal)
# ps = sphere.vertices
# scale_start = Vec3f0[Vec3f0(1,1,rand()) for i=1:length(ps)]
# scale = foldp(scale_gen, scale_start, t)
# colorstart = color_gen(zeros(RGBA{U8}, n), value(t))
# color = foldp(color_gen, colorstart, t)
# rotation = -sphere.normals

# points = visualize((mesh, ps), scale=scale, color=color, rotation=rotation)
cubes = visualize((mesh, position_signal), color=color_signal)

view(cubes, window)

if !isdefined(:runtests)
	renderloop(window)
end
