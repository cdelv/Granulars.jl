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

v = zeros(SVector{3})
w = ones(SVector{3})

println("Unitary v=[0,0,0]")
@btime $unitary($v)

println("Unitary v≢[0,0,0]")
@btime $unitary($w)

println("Check Simulation")
@btime $Check_Simulation($particles)

println("Compute Inertia Tensor")
@btime $Compute_Inertia_Tensor($particles[1])

println("Set Inertia")
@btime $Set_Inertia($particles[1])

println("Angle")
@btime $angle($v, $w)

println("Lab to body")
@btime $Lab_to_body($w,$particles.q[1])

println("Body to lab")
@btime $Body_to_lab($w,$particles.q[1])

