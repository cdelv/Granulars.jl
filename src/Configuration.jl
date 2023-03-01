"""
Struct for representing infinite walls as infinite planes.
- n: Normal vector. It has to be unitary. However, it is normalized in the constructor. 
- Q: A point of the plane. 

- E: Young Modulus.
- G: Shear Modulus. 
- ν: Poisson ratio.

- v: Wall velocity.
- F: Force acting on the wall.
- m: mass of the wall

- static: whether or not to move the wall acording to forces acting on it.
"""
struct Wall
    n::SVector{3, Float64}
    Q::SVector{3, Float64}

    E::Float64
    G::Float64 
    ν::Float64

    v::SVector{3, Float64}
    F::SVector{3, Float64}
    m::Float64

    static::Bool
end

"""
Constructor for Wall. 
For convenience, it normalizes the normal vector automatically.
- n: Normal vector.
- q: A point of the plane. 
- E: Young Modulus.
- ν: Poisson ratio.
"""
function Wall(n::Union{Vector{<:Real}, SVector{3}}, 
    q::Union{Vector{<:Real}, SVector{3}}; 
    E::Real=1.0e6, 
    ν::Real=0.2, 
    v::Union{Vector{<:Real}, SVector{3}}=zeros(3), 
    m::Real=1.0e8,
    static::Bool=true)::Wall
    
    G::Float64 = Float64(E/(2.0*(1.0+ν)))
	N::SVector{3,Float64} = normalize(SVector{3,Float64}(n)) # error if [0,0,0]
	Q::SVector{3,Float64} = SVector{3,Float64}(q)
	Wall(N, Q, Float64(E), G, Float64(ν), SVector{3, Float64}(v), zeros(SVector{3}), Float64(m), static)
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
function Move_r(w::Wall, dt::Float64, cte::Float64=1.0)::Wall
    Set_Q(w, w.Q + w.v*dt*cte)
end

"""
Update velocity according to MD algortihm. 
- v: Velocity vector.
- a: Acceleration vector.
- dt: Time step
- cte: Integration algorithm constant. 
"""
function Move_v(w::Wall, dt::Float64, cte::Float64=1.0)::Wall
    if w.static
        cte = 0.0
    end
    Set_v(w, w.v + w.F*dt*cte/w.m)
end

#=
SET METHODS FOR WALL
=#
function Set_n(w::Wall, n::SVector{3, Float64})::Wall
    @assert norm(n) > 0.0
    return Wall(normalize(n),w.Q,w.E,w.G,w.ν,w.v,w.F,w.m,w.static) # error if [0,0,0]
end
function Set_Q(w::Wall, Q::SVector{3, Float64})::Wall
    return Wall(w.n,Q,w.E,w.G,w.ν,w.v,w.F,w.m,w.static)
end
function Set_E(w::Wall, E::Float64)::Wall
    G::Float64 = E/(2.0*(1.0+w.ν))
    return Wall(w.n,w.Q,E,G,w.ν,w.v,w.F,w.m,w.static)
end
function Set_ν(w::Wall, ν::Float64)::Wall
    G::Float64 = w.E/(2.0*(1.0+ν))
    return Wall(w.n,w.Q,w.E,G,ν,w.v,w.F,w.m,w.static)
end
function Set_v(w::Wall, v::SVector{3, Float64})::Wall
    return Wall(w.n,w.Q,w.E,w.G,w.ν,v,w.F,w.m,w.static)
end
function Set_F(w::Wall, F::SVector{3, Float64})::Wall
    return Wall(w.n,w.Q,w.E,w.G,w.ν,w.v,F,w.m,w.static)
end
function Set_m(w::Wall, m::Float64)::Wall
    @assert m > 0.0
    return Wall(w.n,w.Q,w.E,w.G,w.ν,w.v,w.F,w.m,w.static)
end
function Set_static(w::Wall, static::Bool)::Wall
    return Wall(w.n,w.Q,w.E,w.G,w.ν,w.v,w.F,w.m,static)
end

"""
Struct to configure the simulation. It stores information about the 
time, force parameters, and walls.  
- tf: Maximum simulation time.
- dt: Time step.
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 

- mu: Friction coefficient. TO DO: (Support for multiple species)
- en: Coeficient of restitution. TO DO: (Support for multiple species)

- v: Critical velocity for Thorsten Damping Model.
- thorsten_damping: Whether or not to use Thorsten's Damping Model.

- beam_damping: Whether or not to use damping on beams.
- ζ: Damping ratio for beams. TO DO: (Support for multiple species)
"""
mutable struct Config
	tf::Float64
	dt::Float64
    g::SVector{3, Float64}
    walls::Vector{Wall}

    mu::Float64
    en::Float64

    v::Float64
    thorsten_damping::Bool

    beam_damping::Bool
    ζ::Float64
end

"""
Convenience constructor for Config
- tf: Simulation time.
- dt: Time step.
- g: Gravity. 
- walls: Array of Wall structs that define the simulation bounds. 
- mu: Friction coefficient. TO DO: (Support for multiple species)
- en: Coeficient of restitution. TO DO: (Support for multiple species)
- v: critical velocity for Thorsten Damping Model.
- thorsten_damping: Whether or not to use Thorsten's Damping Model.
- beam_damping: Whether or not to use damping on beams.
- ζ: Damping ratio for beams. TO DO: (Support for multiple species)

Does allocations!!!
"""
function Config(tf::Real, dt::Real; walls::Vector{Wall}=Wall[], 
    g::Union{Vector{<:Real}, SVector{3}}=[0.0,-9.8,0.0], mu::Real=0.4, en::Real=0.9, v::Real=1.0,
    beam_damping::Bool=false, thorsten_damping::Bool=false, ζ::Real=0.05)::Config

	Config(Float64(tf), Float64(dt), SVector{3,Float64}(g), 
        unique(walls), Float64(mu), Float64(en), Float64(v), thorsten_damping, beam_damping, Float64(ζ))
end