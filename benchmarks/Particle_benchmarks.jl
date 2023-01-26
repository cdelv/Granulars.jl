include("../src/Granulars.jl")
using BenchmarkTools

println("CONSTRUCTORS")
println("Default Constructor")
@btime $Particle($zeros(SVector{3}),$zeros(SVector{3}),$zeros(SVector{3}),$Quaternion(1.0I),$zeros(SVector{3}),$zeros(SVector{3}),$0.0,$zeros(SVector{3}),$0.0)
println("First Constructor")
@btime $Particle($[0.0,0.0,0.0], $[0.0,0.0,0.0], $[0.0,0.0,0.0], $1.0)
println("Second Constructor")
@btime $Particle($[0.0,0.0,0.0], $[0.0,0.0,0.0])
println("Third Constructor")
@btime $Particle()
println("")

p = Particle()
dt = 0.1
cte = 0.1

println("MOVE METHODS")
println("Move_r")
@btime $Move_r($p.r,$p.v,$dt,$cte)

println("Move_v")
@btime $Move_v($p.v,$p.a,$dt,$cte)

println("Move_w")
@btime $Move_w($p.w,$p.τ,$p.I,$dt,$cte)

println("Move_q")
@btime $Move_q($p.q,$p.w,$dt)

println("")

println("SET METHODS")
println("Set_r")
@btime $Set_r($p, $zeros(SVector{3}))
println("Set_v")
@btime $Set_v($p, $zeros(SVector{3}))
println("Set_a")
@btime $Set_a($p, $zeros(SVector{3}))
println("Set_q")
@btime $Set_q($p, Quaternion(1.0I))
println("Set_w")
@btime $Set_w($p, $zeros(SVector{3}))
println("Set_τ")
@btime $Set_τ($p, $zeros(SVector{3}))
println("Set_m")
@btime $Set_m($p, $1.0)
println("Set_I")
@btime $Set_I($p, $ones(SVector{3}))
println("Set_rad")
@btime $Set_rad($p, $1.0)
println("")