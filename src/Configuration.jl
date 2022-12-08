"""
Struct for representing infinite walls as planes. They are dimensionally general. 
- D: Dimension.
- n: Normal vector. It has to be unitary. However, it is normalized in the constructor. 
- Q: A point of the plane. 
"""
struct Wall{D}
    n::SVector{D, Float64}
    Q::SVector{D, Float64}
end

"""
Constructor for Wall. For convenience, it normalizes the normal vector automatically.
"""
function Wall(n::Vector{<:Real}, q::Vector{<:Real})
	
	D::Int = length(n)
	
	if D!=length(q)
	    throw(DomainError(q,"n is not the same size as q."))
	end

	N = convert(SVector{D,Float64}, n)
	N = normalize(N)
	Q = convert(SVector{D,Float64}, q)
	Wall(N,Q)
end


"""
Struct to configure the simulation. It stores information about the 
time, force parameters, and walls.  
- D: Dimension.
- tf: Maximum simulation time
- dt: Time step.
- K: Constant for Hertz force calculation. Has to be change to be more general. TO DO: (Support for multiple species)
- gamma: Damping constant
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
- K_kundall: Kundall spring constant.
- mu: Friction coefficient.
"""
struct Config{D}
	tf::Float64
	dt::Float64
	K::Float64
	gamma::Float64
    g::SVector{D, Float64}
    walls::Vector{Wall{D}}
    K_kundall::Float64
    mu::Float64
end

"""
Convenience constructor for Config. It will error if g and Wall have different sizes. 
TO DO: make possible to use no walls.
"""
function Config(tf::Real, dt::Real, g::T; k=4.0e2::Real, gamma=0.2::Real, KK=900.0::Real, mu=0.2::Real, walls=[]::Vector{Wall}) where {T}

	D::Int = length(g)

	TF = convert(Float64,tf)
	DT = convert(Float64,dt)
	K = convert(Float64,k)
	kk = convert(Float64,KK)
	MU = convert(Float64,mu)
	GAMMA = convert(Float64,gamma)
	G = convert(SVector{D,Float64}, g)
	W = convert(Vector{Wall{D}}, unique(walls))


	Config(TF,DT,K,GAMMA,G,W,kk,MU)
end