using CGP
using E4L

CGP.Config.init("cfg/cgp.yml")
@Logging.configure(level=DEBUG)

include("test/individual.jl")
include("test/environment.jl")
include("test/controller.jl")
