using Test
include("../src/Granulars.jl")

# CONSTRUCTORS
q = Particle([0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0], 1.0, rad = 1.0)
Particle([0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0], 1.0)

Particle([0.0,0.0,0.0], [0.0,0.0,0.0])
Particle([0.0,0.0,0.0], [0.0,0.0,0.0], m=1.0)
Particle([0.0,0.0,0.0], [0.0,0.0,0.0], rad=1.0)
Particle([0.0,0.0,0.0], [0.0,0.0,0.0], m=1.0, rad=1.0)

Particle(I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0],I=[0.0,0.0,0.0])
Particle(w=[0.0,0.0,0.0],I=[0.0,0.0,0.0])
Particle(m=1.0,I=[0.0,0.0,0.0])
Particle(rad=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0],I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], m=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], rad=1.0,I=[0.0,0.0,0.0])
Particle(w=[0.0,0.0,0.0], m=1.0,I=[0.0,0.0,0.0])
Particle(w=[0.0,0.0,0.0], rad=1.0,I=[0.0,0.0,0.0])
Particle(m=1.0, rad=1.0)
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], m=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=1.0,I=[0.0,0.0,0.0])
Particle(w=[0.0,0.0,0.0], m=1.0, rad=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], m=1.0, rad=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0],I=[0.0,0.0,0.0])
Particle(w=[0.0,0.0,0.0],I=[0.0,0.0,0.0])
Particle(m=1.0,I=[0.0,0.0,0.0])
Particle(rad=1.0,I=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0])
Particle(r=[0.0,0.0,0.0], m=1.0)
Particle(r=[0.0,0.0,0.0], rad=1.0)
Particle(w=[0.0,0.0,0.0], m=1.0)
Particle(w=[0.0,0.0,0.0], rad=1.0)
Particle(m=1.0, rad=1.0)
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], m=1.0)
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=1.0)
Particle(w=[0.0,0.0,0.0], m=1.0, rad=1.0)
Particle(r=[0.0,0.0,0.0], w=[0.0,0.0,0.0], m=1.0, rad=1.0)

p = Particle(r=[1.0,1.0,1.0],v=[1.0,1.0,1.0],w=[1.0,1.0,1.0])
p = Set_I(p,SVector(1.0,1.0,1.0))
dt::Float64 = 1.0
cte::Float64 = 1.0

@test p.r == SVector{3,Float64}([1,1,1])
#@test π ≈ 3.14 atol=0.01

# MOVE METHODS
@test SVector{3,Float64}([2,2,2]) == Move_r(p.r,p.v,dt,cte)
@test SVector{3,Float64}([2,2,2]) == Move_v(p.r,p.v,dt,cte)
@test SVector{3,Float64}([2,2,2]) == Move_w(p.r,p.v,p.I,dt,cte) # test for anisotropic I
#@test Move_q(p.q,p.w,dt)) how to test this? w=0?

# SET METHODS
p = Particle()
@test Set_r(p, SVector{3,Float64}([2,2,2])).r == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_v(p, SVector{3,Float64}([2,2,2])).v == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_a(p, SVector{3,Float64}([2,2,2])).a == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_q(p, Quaternion(2.0I)).q == Quaternion(2.0I)
p = Particle()
@test Set_w(p, SVector{3,Float64}([2,2,2])).w == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_τ(p, SVector{3,Float64}([2,2,2])).τ == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_m(p, 2.0).m == 2.0
p = Particle()
@test Set_I(p, SVector{3,Float64}([2,2,2])).I == SVector{3,Float64}([2,2,2])
p = Particle()
@test Set_rad(p, 2.0).rad == 2.0