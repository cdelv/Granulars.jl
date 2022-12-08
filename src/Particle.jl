# @forward Citizen.person name --> Lazy.@forward  for inheritance
"""
Abstract type for particles and inheritance. All particles should have: 
- D: Dimension.
- r: Position vector.
- v: Velocity vector.
- a: Acceleration vector.
- m: Mass. In the future, it may change for density. 
"""
abstract type AbstractParticle{D} end

"""
Struct to represent particles of dimension D. For example, 2D or 3D.
It's immutable for performance. Inspired from AtomsBase.jl. It uses acceleration, that
way, one can avoid dividing by the mass every step.
- r is the position of the center of mass of the particle.
- v is the velocity of the center of mass of the particle.
- a is the acceleration of the center of mass of the particle.
- m is the mass of the particle.
- rad is the radius of the particle.
"""
struct Particle{D} <: AbstractParticle{D}
    r::SVector{D, Float64}
    v::SVector{D, Float64}
    a::SVector{D, Float64}
    m::Float64
    rad::Float64
end

#=
 CONSTRUCTORS -> MAY BE IMPROVE?
=#
"""
These constructors are just for convenience. This way, one doesn't need to specify 
every parameter to create a particle. However, they aren't type-stable. 

Example -> Particle([0,0,0],[0,0,0],1,rad=1)
"""
function Particle(r::Vector{<:Real}, v::Vector{<:Real}, m::Real; rad::Real=1.0)::Particle
    
    n = length(r)

    if n!=length(v)
        throw(DomainError(v,"v is not the same size as r"))
    end

    Particle(SVector{n,Float64}(r),SVector{n,Float64}(v),zeros(SVector{n}),Float64(m),Float64(rad))
end

"""
Example -> Particle([0,0,0],[0,0,0],m=1,rad=1)
"""
function Particle(r::Vector{<:Real}, v::Vector{<:Real}; m::Real=1.0, rad::Real=1.0)::Particle
    
    n = length(r)
    
    if n!=length(v)
        throw(DomainError(v,"v is not the same size as r"))
    end
    
    Particle(SVector{n,Float64}(r),SVector{n,Float64}(v),zeros(SVector{n}),Float64(m),Float64(rad))
end

"""
Example -> Particle(r=[0,0,0],m=1)
Example -> Particle().
"""
function Particle(;r::Vector{<:Real}=[0.0,0.0,0.0], m::Real=1.0, rad::Real=1.0)::Particle

    n = length(r)

    Particle(SVector{n,Float64}(r),SVector{n,Float64}(zeros(n)),SVector{n,Float64}(zeros(n)),Float64(m),Float64(rad))
end

#=
 METHODS
=#
"""
Get the dimension of a particle P.
"""
Get_Dim(P::AbstractParticle{D}) where {D} = D

"""
For Newton's equations of motion integration, we use PEFRL. Most methods like
Leapfrog, Verlet, Forest-Ruth, etc. work similarly. These methods allow 
changing the integration algorithm quickly.
- r: Position vector.
- v: Velocity vector.
- dt: Time step
- cte: Integration algorithm constant. 
"""
function Move_r(r::SVector{D, Float64}, v::SVector{D, Float64}, dt::Float64, cte=1.0::Float64)::SVector{D, Float64} where {D}
    return r + v*dt*cte
end

"""
Update velocity according to MD algortihm. 
- v: Velocity vector.
- a: Acceleration vector.
- dt: Time step
- cte: Integration algorithm constant. 
"""
function Move_v(v::SVector{D, Float64}, a::SVector{D, Float64}, dt::Float64, cte=1.0::Float64)::SVector{D, Float64} where {D}
    return v + a*dt*cte
end