using Test
include("../src/Granulars.jl")

v = zeros(SVector{3})
w = ones(SVector{3})

@test unitary(v) == v 
@test unitary(w) == w/sqrt(3)