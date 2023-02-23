# @forward Citizen.person name --> Lazy.@forward  for inheritance
"""
Abstract type for particles and inheritance. All particles should have: 
- r: position of the center of mass of the particle.
- v: velocity of the center of mass of the particle.
- a: acceleration of the center of mass of the particle.
- q: quaternion that represents the orientation of the particle .
- w: angular velocity of the particle in the body frame.
- τ: torque over the particle in the body frame.
- m: is the mass of the particle. Change for density?
- I: inertia of the particle in the principal axes reference frame. 
- rad: radius of the particle.
"""
abstract type AbstractParticle end

"""
Struct to represent particles. It's immutable for performance. Inspired from AtomsBase.jl. 
I use acceleration, that way, I can avoid dividing by the mass every step. 
Due to the rotation algorithm, I use torque as the inertia is used in a different way than mass. 
- r: position of the center of mass of the particle.
- v: velocity of the center of mass of the particle.
- a: acceleration of the center of mass of the particle.
- q: quaternion that represents the orientation of the particle.
- w: angular velocity of the particle in de body frame.
- τ: torque over the particle in de body frame.
- m: is the mass of the particle. Change for density?
- I: inertia of the particle in the principal axes reference frame. 
- rad: radius of the particle.
- E: Young Modulus.
- G: Shear Modulus. 
- ν: Poisson ratio.
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

    E::Float64
    G::Float64 
    ν::Float64
end

#=
 CONSTRUCTORS
=#
"""
These constructors are just for convenience. This way, one doesn't need to specify 
every parameter to create a particle.

Example -> Particle([0,0,0],[0,0,0],1,rad=1)
"""
function Particle(r::Union{Vector{<:Real}, SVector{3}}, 
    v::Union{Vector{<:Real}, SVector{3}}, 
    w::Union{Vector{<:Real}, SVector{3}}, 
    m::Real; 
    rad::Real=1.0, 
    E::Real=1.0e6, 
    ν::Real=0.2)::Particle

    G::Float64 = Float64(E/(2.0*(1.0+ν)))
    p::Particle = Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        SVector{3,Float64}(w),
        zeros(SVector{3}),
        Float64(m),
        ones(SVector{3}),
        Float64(rad),
        Float64(E),
        G, Float64(ν))
    p = Set_Inertia(p) # Adds the innertia tensor and particle orientation. Defined in Utils.jl
    Set_w(p,Lab_to_body(w,p.q))
end

"""
Example -> Particle([0,0,0],[0,0,0],m=1,rad=1)
"""
function Particle(r::Union{Vector{<:Real}, SVector{3}}, 
    v::Union{Vector{<:Real}, SVector{3}}; 
    m::Real=1.0, 
    rad::Real=1.0, 
    E::Real=1.0e6, 
    ν::Real=0.2)::Particle

    G::Float64 = Float64(E/(2.0*(1.0+ν)))
    p::Particle = Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        zeros(SVector{3}),
        zeros(SVector{3}),
        Float64(m),
        ones(SVector{3}),
        Float64(rad),
        Float64(E),
        G, Float64(ν))
    Set_Inertia(p) # Adds the innertia tensor and particle orientation. Defined in Utils.jl
end

"""
This one has alocations and its slow

Example -> Particle(r=[0,0,0],m=1)
Example -> Particle().
"""
function Particle(;r::Union{Vector{<:Real}, SVector{3}}=[0.0,0.0,0.0], 
    v::Union{Vector{<:Real}, SVector{3}}=[0.0,0.0,0.0], 
    w::Union{Vector{<:Real}, SVector{3}}=[0.0,0.0,0.0], 
    m::Real=1.0, 
    I::Union{Vector{<:Real}, SVector{3}}=[1.0,1.0,1.0], 
    rad::Real=1.0, 
    E::Real=1.0e6, 
    ν::Real=0.2)::Particle

    G::Float64 = Float64(E/(2.0*(1.0+ν)))
    p::Particle = Particle(SVector{3,Float64}(r),
        SVector{3,Float64}(v),
        zeros(SVector{3}),
        Quaternion(1.0I),
        SVector{3,Float64}(w),
        zeros(SVector{3}),
        Float64(m),
        SVector{3,Float64}(I),
        Float64(rad),
        Float64(E),
        G, Float64(ν))
    p = Set_Inertia(p) # Adds the innertia tensor and particle orientation. Defined in Utils.jl
    Set_w(p,Lab_to_body(p.w,p.q))
end

#=
METHODS FOR MD
=# 
"""
For Newton's equations of motion integration.
Most methods like Leapfrog, Verlet, Forest-Ruth, etc. work similarly. 
These method allows changing the integration algorithm quickly.
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
function Move_w(w::SVector{3, Float64}, τ::SVector{3, Float64}, II::SVector{3, Float64}, dt::Float64, cte=1.0::Float64, niter=3::Int64)::SVector{3, Float64}
    # Initial guess as recomended in the article
    ww::SVector{3, Float64} = w
    
    # Solve non-linear system
    for i in 1:niter
        @inbounds ww = w + SVector(
            τ[1] + 0.5*(w[2]*w[3]+ww[2]*ww[3])*(II[2]-II[3]),
            τ[2] + 0.5*(w[3]*w[1]+ww[3]*ww[1])*(II[3]-II[1]),
            τ[3] + 0.5*(w[1]*w[2]+ww[1]*ww[2])*(II[1]-II[2]))*dt*cte./II
    end

    return ww
end

"""
Update orientation according to MD algortihm. I use: 

- Algorithm for numerical integration of the rigid-body equations of motion, Igor P. Omelyan, 1998

- q: quaternion that represents the orientation.
- w: Angular velocity vector.
- dt: Time step.
- cte: Integration algorithm constant. 
"""
function Move_q(q::Quaternion{Float64}, w::SVector{3, Float64}, dt::Float64)::Quaternion{Float64}
    Q_q_dt::Quaternion{Float64} = dquat(q, w)
    a1::Float64 = 1.0 - dt*dt*dot(w,w)/16.0
    a2::Float64 = 1.0 + dt*dt*dot(w,w)/16.0
    return (a1*q + Q_q_dt*dt)/a2
end

#=
SET METHODS FOR PARTICLE
=#
function Set_r(p::Particle, r::SVector{3, Float64})::Particle
    return Particle(r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_v(p::Particle, v::SVector{3, Float64})::Particle
    return Particle(p.r,v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_a(p::Particle, a::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,a,p.q,p.w,p.τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_q(p::Particle, q::Quaternion{Float64})::Particle
    return Particle(p.r,p.v,p.a,q,p.w,p.τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_w(p::Particle, w::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,p.a,p.q,w,p.τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_τ(p::Particle, τ::SVector{3, Float64})::Particle
    return Particle(p.r,p.v,p.a,p.q,p.w,τ,p.m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_m(p::Particle, m::Float64)::Particle
    @assert m > 0.0
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,m,p.I,p.rad,p.E,p.G,p.ν)
end
function Set_I(p::Particle, I::SVector{3, Float64})::Particle
    @assert all(>(0.0), I)
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,I,p.rad,p.E,p.G,p.ν)
end
function Set_rad(p::Particle, rad::Float64)::Particle
    @assert rad > 0.0
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,rad,p.E,p.G,p.ν)
end
function Set_E(p::Particle, E::Float64)::Particle
    G::Float64 = E/(2.0*(1.0+ν))
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad,E,G,p.ν)
end
function Set_ν(p::Particle, ν::Float64)::Particle
    G::Float64 = E/(2.0*(1.0+ν))
    return Particle(p.r,p.v,p.a,p.q,p.w,p.τ,p.m,p.I,p.rad,p.E,G,ν)
end