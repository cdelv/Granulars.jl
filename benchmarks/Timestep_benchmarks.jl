include("../src/Granulars.jl")
using BenchmarkTools

W1 = Wall([1,0,0],[0,0,0])
W2 = Wall([-1,0,0],[25,0,0])
W3 = Wall([0,1,0],[0,0,0])
W4 = Wall([0,-1,0],[0,25,0])
W5 = Wall([0,0,1],[0,0,0])
W6 = Wall([0,0,-1],[0,0,25])
walls = [W1,W2,W3,W4,W5,W6]
C = Config(10.0,0.1,walls=walls)
data = [Particle((25-2.0)*rand(3)+[1,1,1],4.0*rand(3)) for i = 1:100]
particles = StructArray(data)
Cutoff::Float64 = 4*maximum(particles.rad)
system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
list = neighborlist!(system) # Type Warning
kundall_particles = ExtendableSparseMatrix(zeros(length(data), length(data)))
kundall_walls = ExtendableSparseMatrix(zeros(length(data), length(C.walls)))

println("time step")
@btime time_step($particles, $C, $list, $kundall_particles, $kundall_walls)