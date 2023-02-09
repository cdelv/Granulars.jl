"""
Struct for representing infinite walls as planes. They are dimensionally general. 
- D: Dimension.
- n: Normal vector. It has to be unitary. However, it is normalized in the constructor. 
- Q: A point of the plane. 
"""
struct Wall
    n::SVector{3, Float64}
    Q::SVector{3, Float64}
end

"""
Constructor for Wall. 
For convenience, it normalizes the normal vector automatically.
- n: Normal vector.
- q: A point of the plane. 
"""
function Wall(n::Vector{<:Real}, q::Vector{<:Real})::Wall
	N = convert(SVector{3,Float64}, n)
	N = normalize(N) # error if [0,0,0]
	Q = convert(SVector{3,Float64}, q)
	Wall(N,Q)
end


"""
Struct to configure the simulation. It stores information about the 
time, force parameters, and walls.  
- tf: Maximum simulation time.
- dt: Time step.
- K: Constant for Hertz force calculation. CHANGE FOR EXACT VAL. TO DO: (Support for multiple species)
- gamma: Damping constant.
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- K_kundall: Kundall spring constant. TO DO: (Support for multiple species)
- mu: Friction coefficient. TO DO: (Support for multiple species)
- E: Young Modulus.
- G: Shear Modulus. 
"""
struct Config
	tf::Float64
	dt::Float64
	K::Float64
	gamma::Float64
    g::SVector{3, Float64}
    walls::Vector{Wall}
    K_cundall::Float64
    mu::Float64
    G::Float64
    E::Float64
end

"""
Convenience constructor for Config. It will error if g and Wall have different sizes. 
- tf: Maximum simulation time
- dt: Time step.
- K: Constant for Hertz force calculation. Has to be change to be more general.
- gamma: Damping constant
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- K_kundall: Kundall spring constant.
- mu: Friction coefficient.
- E: Young Modulus.
- G: Shear Modulus. 

Does allocations!!!
"""
function Config(tf::Real, dt::Real; walls::Vector{Wall}=Wall[], 
    g::Vector{<:Real}=[0,-9.8,0], k::Real=5.0e3, 
    gamma::Real=500, KK::Real=900.0, mu::Real=0.4,
    G::Real=120000.0, E::Real=80000.0)::Config

	Config(convert(Float64,tf), convert(Float64,dt), convert(Float64,k),
        convert(Float64,gamma), convert(SVector{3,Float64}, g), 
        convert(Vector{Wall}, unique(walls)),convert(Float64,KK), convert(Float64,mu),
        convert(Float64,G), convert(Float64,E))
end