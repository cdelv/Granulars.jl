"""
Struct for representing infinite walls as planes. They are dimensionally general. 
- D: Dimension.
- n: Normal vector. It has to be unitary. However, it is normalized in the constructor. 
- Q: A point of the plane. 

- E: Young Modulus.
- G: Shear Modulus. 
- ν: Poisson ratio
"""
struct Wall
    n::SVector{3, Float64}
    Q::SVector{3, Float64}

    E::Float64
    G::Float64 
    ν::Float64
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
    G::Real=1.0e6)::Wall
    
    ν::Float64 = Float64(E/(2.0*G) - 1.0) # Poisson ratio
	N::SVector{3,Float64} = normalize(SVector{3,Float64}(n)) # error if [0,0,0]
	Q::SVector{3,Float64} = SVector{3,Float64}(q)
	Wall(N,Q,Float64(E),Float64(G),ν)
end

#=
SET METHODS
=#
function Set_n(w::Wall, n::SVector{3, Float64})::Wall
    return Wall(normalize(n),w.Q,w.E,w.G,w.ν)
end
function Set_Q(w::Wall, Q::SVector{3, Float64})::Wall
    return Wall(w.n,Q,w.E,w.G,w.ν)
end
function Set_E(w::Wall, E::Float64)::Wall
    ν::Float64 = E/(2.0*w.G) - 1.0
    return Wall(w.n,w.Q,E,w.G,ν)
end
function Set_G(w::Wall, G::Float64)::Wall
    ν::Float64 = w.E/(2.0*G) - 1.0
    return Wall(w.n,w.Q,w.E,G,ν)
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
"""
struct Config
	tf::Float64
	dt::Float64
    g::SVector{3, Float64}
    walls::Vector{Wall}

    mu::Float64
    en::Float64
end

"""
Convenience constructor for Config
- tf: Maximum simulation time
- dt: Time step.
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- mu: Friction coefficient.

Does allocations!!!
"""
function Config(tf::Real, dt::Real; walls::Vector{Wall}=Wall[], 
    g::Vector{<:Real}=[0.0,-9.8,0.0], mu::Real=0.4, en::Real=0.9)::Config

	Config(Float64(tf), Float64(dt), SVector{3,Float64}(g), 
        Vector{Wall}(unique(walls)), Float64(mu), Float64(en))
end