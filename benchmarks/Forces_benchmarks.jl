include("../src/Granulars.jl")
using BenchmarkTools

W1 = Wall([1,0,0],[0,0,0])
W2 = Wall([-1,0,0],[25,0,0])
W3 = Wall([0,1,0],[0,0,0])
W4 = Wall([0,-1,0],[0,25,0])
W5 = Wall([0,0,1],[0,0,0])
W6 = Wall([0,0,-1],[0,0,25])
walls = [W1,W2,W3,W4,W5,W6]
C = Config(10.0,0.1,walls)
data = [Particle((25-2.0)*rand(3)+[1,1,1],4.0*rand(3)) for i = 1:100]
particles = StructArray(data)
Cutoff::Float64 = 4*maximum(particles.rad)
system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
list = neighborlist!(system) # Type Warning
kundall_particles = ExtendableSparseMatrix(zeros(length(data), length(data)))
kundall_walls = ExtendableSparseMatrix(zeros(length(data), length(C.walls)))

println("Kundall_friction")
@btime Kundall_friction($1.0, $100.0, $C)
println("Damping force")
@btime Damping_Force($0.05, $0.5, $12.0, $C)
println("Hertz force")
@btime Hertz_Force($0.05, $C)

println("Force with walls")
@btime Force_With_Walls($particles, $50, $C, $kundall_walls)

println("Force with pairs")
@btime Force_With_Pairs($particles, $C, $list, $kundall_particles)

println("Calculate Forces")
@btime Calculate_Forces($particles, $list, $C, $kundall_particles, $kundall_walls)