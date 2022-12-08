#Usefull macros
# @time @btime @benchmark @allocated @code_warntype @simd @turbo @inbounds

using StaticArrays
using StructArrays
using CellListMap
using CSV
using LinearAlgebra
using Accessors
using SparseArrays: dropzeros!
#using SuiteSparseGraphBLAS # explore this library after the updates that remove elements
using ExtendableSparse: ExtendableSparseMatrix


include("Particle.jl")
include("Configuration.jl")
include("Time_integration.jl")
include("Writter.jl")
include("Forces.jl")


# Documentation
# Check simulation state
# Function to initialize particles
# 

# 3. Implement Rotations
# 4. Implement BeamForces
# 5. Implement MultiSpheres
# 6. Restart simulation and start from file
;