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
Constructor for Wall. For convenience, it normalizes the normal vector automatically.
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
- tf: Maximum simulation time
- dt: Time step.
- K: Constant for Hertz force calculation. Has to be change to be more general. TO DO: (Support for multiple species)
- gamma: Damping constant
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- K_kundall: Kundall spring constant. TO DO: (Support for multiple species)
- mu: Friction coefficient. TO DO: (Support for multiple species)
"""
struct Config
	tf::Float64
	dt::Float64
	K::Float64
	gamma::Float64
    g::SVector{3, Float64}
    walls::Vector{Wall}
    K_kundall::Float64
    mu::Float64
end

"""
Convenience constructor for Config. It will error if g and Wall have different sizes. 
TO DO: make possible to use no walls.
- tf: Maximum simulation time
- dt: Time step.
- K: Constant for Hertz force calculation. Has to be change to be more general.
- gamma: Damping constant
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- K_kundall: Kundall spring constant.
- mu: Friction coefficient.

Does allocations!!!
"""
function Config(tf::Real, dt::Real, walls::Vector{Wall}; g::Vector{<:Real}=[0,-9.8,0], k=1.0e5::Real, gamma=500::Real, KK=500.0::Real, mu=1.2::Real)::Config
	TF = convert(Float64,tf)
	DT = convert(Float64,dt)
	K = convert(Float64,k)
	kk = convert(Float64,KK)
	MU = convert(Float64,mu)
	GAMMA = convert(Float64,gamma)
	G = convert(SVector{3,Float64}, g)
	W = convert(Vector{Wall}, unique(walls)) # remove reapeted walls
	Config(TF,DT,K,GAMMA,G,W,kk,MU)
end