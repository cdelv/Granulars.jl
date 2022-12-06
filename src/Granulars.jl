#Usefull macros
# @time @btime @benchmark @allocated @code_warntype @simd @turbo @inbounds

using StaticArrays
using StructArrays
using CellListMap
using CSV
using LinearAlgebra

include("Particle.jl")
include("Configuration.jl")
include("Time_integration.jl")
include("Writter.jl")
include("Forces.jl")


# 1. Implement Friction and Damping Forces
# 2. Test Everithing - Energy and other metrics (check no particles left the box)
# 3. Implement Rotations
# 4. Implement BeamForces
# 5. Implement MultiSpheres
# 6. Restart simulation and start from file
;