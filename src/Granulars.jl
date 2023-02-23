#Usefull macros
# @time @btime @benchmark @allocated @code_warntype @simd @turbo @inbounds
using StaticArrays
using StructArrays
using CellListMap
using CSV
using LinearAlgebra
using ReferenceFrameRotations
using Roots
#using Accessors
#using SparseArrays: dropzeros!
#using SuiteSparseGraphBLAS # explore this library after the updates that remove elements
#using ExtendableSparse: ExtendableSparseMatrix

include("Particle.jl")
include("Configuration.jl")
include("Beams.jl")
include("Time_integration.jl")
include("Writter.jl")
include("Forces.jl")
include("Utils.jl")

"""
Tests Ideas:
    - Use ν instead of G
    - Re run all examples
    - Complete documentation
    - Revisit friction force and shear tangent force
    - Thorsten Damping model
    - Fracture Model 
    - Time step stimate
"""


# Delete old data files
# Check if directory exists

# 6. Restart simulation and start from file

# 8. Calculate quantities like energy.
# 9. Create Docs page.
;