#Usefull macros
# @time @btime @benchmark @allocated @code_warntype @simd @turbo @inbounds
using StaticArrays
using StructArrays
using CellListMap
using CSV
using LinearAlgebra
#using Accessors
using SparseArrays: dropzeros!
#using SuiteSparseGraphBLAS # explore this library after the updates that remove elements
using ExtendableSparse: ExtendableSparseMatrix
using ReferenceFrameRotations
using Distributions

include("Particle.jl")
include("Configuration.jl")
include("Beams.jl")
include("Time_integration.jl")
include("Writter.jl")
include("Forces.jl")
include("Utils.jl")

"""
Tests Ideas:
    - Maxwell Boltzman (ideal gas 3D)
    - No Particles outside the box
    - Cambiar constantes de Hertz y damping por lo que es
    - Normalized quaternions
    - viga cantilever
    - viga de flexion
    - Oscilación viga empotrada
    - sol analítica trompo
    - Agradecimiento a ricardo
    - rotar esferas y ver fuerza 0
    - T que da vueltas
"""

# Function to fix particles

# Delete old data files
# Check if directory exists

# 6. Restart simulation and start from file

# 8. Calculate quantities like energy.
# 9. Create Docs page.
;