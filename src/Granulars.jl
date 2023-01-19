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
include("Time_integration.jl")
include("Writter.jl")
include("Forces.jl")
include("Utils.jl")

"""
Tests Ideas:
    - Maxwell Boltzman (ideal gas 3D)
    - No Particles outside the box
    - Cambiar constantes de Hertz y damping por lo que es
"""

"""
DUDAS:
1. Como inicializar la orientación de los ejes principales?
2. 
"""
# Cundall not Kundall, FIX

# Delete old data files
# Check if directory exists

# 3. Compute Inertia tensor and principal axys
# 
# 4. Implement BeamForces
# 5. Implement MultiSpheres and sphere fitting
# - https://jekel.me/2015/Least-Squares-Sphere-Fit/
# - https://stackoverflow.com/questions/70187153/how-can-i-fit-a-sphere-to-a-set-of-3d-points-with-least-squares-method
# - https://github.com/cserteGT3/RANSAC.jl

# 6. Restart simulation and start from file

# 8. Calculate quantities like energy.
# 9. Create Docs page.
;