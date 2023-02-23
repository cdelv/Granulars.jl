"""
Struct for representing infinite walls as infinite planes.
- n: Normal vector. It has to be unitary. However, it is normalized in the constructor. 
- Q: A point of the plane. 

- E: Young Modulus.
- G: Shear Modulus. 
- ν: Poisson ratio
- F: Force acting on the wall
"""
struct Wall
    n::SVector{3, Float64}
    Q::SVector{3, Float64}

    E::Float64
    G::Float64 
    ν::Float64

    F::SVector{3, Float64}
end

"""
Constructor for Wall. 
For convenience, it normalizes the normal vector automatically.
- n: Normal vector.
- q: A point of the plane. 
"""
function Wall(n::Union{Vector{<:Real}, SVector{3}}, 
    q::Union{Vector{<:Real}, SVector{3}}; 
    E::Real=1.0e6, 
    ν::Real=0.2)::Wall
    
    G::Float64 = Float64(E/(2.0*(1.0+ν)))
	N::SVector{3,Float64} = normalize(SVector{3,Float64}(n)) # error if [0,0,0]
	Q::SVector{3,Float64} = SVector{3,Float64}(q)
	Wall(N,Q,Float64(E),G,Float64(ν),zeros(SVector{3}))
end

#=
SET METHODS FOR WALL
=#
function Set_n(w::Wall, n::SVector{3, Float64})::Wall
    return Wall(normalize(n),w.Q,w.E,w.G,w.ν,w.F) # error if [0,0,0]
end
function Set_Q(w::Wall, Q::SVector{3, Float64})::Wall
    return Wall(w.n,Q,w.E,w.G,w.ν,w.F)
end
function Set_E(w::Wall, E::Float64)::Wall
    G::Float64 = E/(2.0*(1.0+w.ν))
    return Wall(w.n,w.Q,E,G,w.ν,w.F)
end
function Set_ν(w::Wall, ν::Float64)::Wall
    G::Float64 = w.E/(2.0*(1.0+ν))
    return Wall(w.n,w.Q,w.E,G,ν,w.F)
end
function Set_F(w::Wall, F::SVector{3, Float64})::Wall
    return Wall(w.n,w.Q,w.E,w.G,w.ν,F)
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

- v: critical velocity for Thorsten Damping Model.
"""
mutable struct Config
	tf::Float64
	dt::Float64
    g::SVector{3, Float64}
    walls::Vector{Wall}

    mu::Float64
    en::Float64

    v::Float64
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

Does allocations!!!
"""
function Config(tf::Real, dt::Real; walls::Vector{Wall}=Wall[], 
    g::Union{Vector{<:Real}, SVector{3}}=[0.0,-9.8,0.0], mu::Real=0.4, en::Real=0.9, v::Real=1.0)::Config

	Config(Float64(tf), Float64(dt), SVector{3,Float64}(g), 
        unique(walls), Float64(mu), Float64(en), Float64(v))
end