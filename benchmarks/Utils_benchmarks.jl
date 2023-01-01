include("../src/Granulars.jl")
using BenchmarkTools

v = zeros(SVector{3})
w = ones(SVector{3})

println("Unitary v=[0,0,0]")
@btime unitary($v)
println("Unitary v≢[0,0,0]")
@btime unitary($w)