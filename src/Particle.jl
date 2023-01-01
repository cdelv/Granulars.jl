# @forward Citizen.person name --> Lazy.@forward  for inheritance
"""
Abstract type for particles and inheritance. All particles should have: 
- r: position of the center of mass of the particle.
- v: velocity of the center of mass of the particle.
- a: acceleration of the center of mass of the particle.
- q: quaternion that represents the orientation of the particle 
- w: angular velodity of the particle 
- τ: torque over the particle
- m: is the mass of the particle. Change for density?
- I: inertia of the particle in the principal axes reference frame. 
- rad: radius of the particle.
"""
abstract type AbstractParticle end

"""
Struct to represent particles. It's immutable for performance. Inspired from AtomsBase.jl. 
I use acceleration, that way, one can avoid dividing by the mass every step. Due to the rotation algorithm,
I use torque as the inertia is used in a different way. 
- r: position of the center of mass of the particle.
- v: velocity of the center of mass of the particle.
- a: acceleration of the center of mass of the particle.
- q: quaternion that represents the orientation of the particle 
- w: angular velodity of the particle 
- τ: torque over the particle
- m: is the mass of the particle. Change for density?
- I: inertia of the particle in the principal axes reference frame. 
- rad: radius of the particle.
"""
struct Particle <: AbstractParticle
    r::SVector{3, Float64}
    v::SVector{3, Float64}
    a::SVector{3, Float64}

    q::Quaternion{Float64}
    w::SVector{3, Float64}
    τ::SVector{3, Float64}

    m::Float64
    I::SVector{3, Float64}

    rad::Float64
end

#=
 CONSTRUCTORS
=#
"""
These constructors are just for convenience. This way, one doesn't need to specify 
every parameter to create a particle. However, they aren't type-stable. 

Example -> Particle([0,0,0],[0,0,0],1,rad=1)
"""
function Particle(r::Vector{<:Real}, v::Vector{<:Real}, w::Vector{<:Real}, m::Real; rad::Real=1.0)::Particle
    Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        SVector{3,Float64}(w),
        zeros(SVector{3}),
        Float64(m),
        ones(SVector{3}),
        Float64(rad))
end

"""
Example -> Particle([0,0,0],[0,0,0],m=1,rad=1)
"""
function Particle(r::Vector{<:Real}, v::Vector{<:Real}; m::Real=1.0, rad::Real=1.0)::Particle
    Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        zeros(SVector{3}),
        zeros(SVector{3}),
        Float64(m),
        ones(SVector{3}),
        Float64(rad))
end

"""
This one has alocations and its slow

Example -> Particle(r=[0,0,0],m=1)
Example -> Particle().
"""
function Particle(;r::Vector{<:Real}=[0.0,0.0,0.0], v::Vector{<:Real}=[0.0,0.0,0.0], w::Vector{<:Real}=[0.0,0.0,0.0], m::Real=1.0, I::Vector{<:Real}=[1.0,1.0,1.0], rad::Real=1.0)::Particle
    Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        SVector{3,Float64}(w),
        zeros(SVector{3}),
        Float64(m),
        SVector{3,Float64}(I),
        Float64(rad))
end

#=
METHODS FOR MD
=# 
"""
For Newton's equations of motion integration, I use PEFRL:

- Optimized Forest–Ruth- and Suzuki-like algorithms for integration of motion in many-body systems, I.P. Omelyan, I.M. Mryglodab and R. Folk, 2002

Most methods like Leapfrog, Verlet, Forest-Ruth, etc. work similarly. These methods allow 
changing the integration algorithm quickly.
- r: Position vector.
- v: Velocity vector.
- dt: Time step.
- cte: Integration algorithm constant. 
"""
function Move_r(r::SVector{3, Float64}, v::SVector{3, Float64}, dt::Float64, cte=1.0::Float64)::SVector{3, Float64}
    return r + v*dt*cte
end

"""
Update velocity according to MD algortihm. 
- v: Velocity vector.
- a: Acceleration vector.
- dt: Time step
- cte: Integration algorithm constant. 
"""
function Move_v(v::SVector{3, Float64}, a::SVector{3, Float64}, dt::Float64, cte=1.0::Float64)::SVector{3, Float64}
    return v + a*dt*cte
end

"""
Update angular velocity according to MD algortihm. I use: 

- Algorithm for numerical integration of the rigid-body equations of motion, Igor P. Omelyan, 1998

- w: Angular velocity vector.
- τ: Torque vector.
- dt: Time step
- cte: Integration algorithm constant. 
"""
function Move_w(w::SVector{3, Float64}, τ::SVector{3, Float64}, II::SVector{3, Float64}, dt::Float64, cte=1.0::Float64)::SVector{3, Float64}
    return w + dt*cte*SVector(
        τ[1] + w[2]*w[3]*(II[2]-II[3]),
        τ[2] + w[3]*w[1]*(II[3]-II[1]),
        τ[3] + w[1]*w[2]*(II[1]-II[2]))./II
end

"""
Update orientation according to MD algortihm. 
- q: quaternion that represents the orientation.
- w: Angular velocity vector.
- dt: Time step.
- cte: Integration algorithm constant. 
"""
function Move_q(q::Quaternion{Float64}, w::SVector{3, Float64}, dt::Float64)::Quaternion{Float64}
    #angle_to_quat(0.5, 0, 0, :XYZ) and quat_to_angle(q::Quaternion, :ZYX)
    Q_q_dt::Quaternion{Float64} = dquat(q, w)*dt
    a1::Float64 = 1.0 - dt*dt*dot(w,w)/16
    a2::Float64 = 1.0 + dt*dt*dot(w,w)/16
    return (a1*q + Q_q_dt)/a2
end

#=
SET METHODS
=#
function Set_r(p::Particle, r::SVector{3, Float64})::Particle # Convert svector, acept normal vector
    return Particle(r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad)
end
function Set_v(p::Particle, v::SVector{3, Float64})::Particle
    return Particle(p.r,v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad)
end
function Set_a(p::Particle, a::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,a,p.q,p.w,p.τ,p.m,p.I,p.rad)
end
function Set_q(p::Particle, q::Quaternion{Float64})::Particle
    return Particle(p.r,p.v,p.a,q,p.w,p.τ,p.m,p.I,p.rad)
end
function Set_w(p::Particle, w::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,p.a,p.q,w,p.τ,p.m,p.I,p.rad)
end
function Set_τ(p::Particle, τ::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,p.a,p.q,p.w,τ,p.m,p.I,p.rad)
end
function Set_m(p::Particle, m::Float64)::Particle
    @assert m > 0.0
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,m,p.I,p.rad)
end
function Set_I(p::Particle, I::SVector{3, Float64})::Particle
    @assert all(>(0.0), I)
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,I,p.rad)
end
function Set_rad(p::Particle, rad::Float64)::Particle
    @assert rad > 0.0
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,rad)
end