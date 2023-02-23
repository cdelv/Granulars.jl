include("../src/Granulars.jl")
using BenchmarkTools

println("WALL CONSTRUCTORS")
println("Default constructor")
@btime $Wall($ones(SVector{3}),$ones(SVector{3}))
println("First constructor")
@btime $Wall($[1,1,1],$[1,1,1])
println("")

println("CONFIG CONSTRUCTORS")
W = Wall([1,1,1],[1,1,1])
println("Default constructor")
@btime $Config($10.0,$0.1,$1.0e6,$500.0,$ones(SVector{3}),$[W],$900.0,$0.4)
println("First constructor")
@btime $Config($10.0,$0.1)
println("")