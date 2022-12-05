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
Constructor for Wall. For convenience, it normalizes the normal vector.
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
- K: Constant for Hertz force calculation. Has to be change to be more general. 
(Support for multiple species)
- g: Simulation gravity vector. 
- walls: Array of Wall structs that define the simulation bounds. 
"""
struct Config{D}
	tf::Float64
	dt::Float64
	K::Float64
    g::SVector{D, Float64}
    walls::Vector{Wall{D}}
end

"""
Convenience constructor for Config. It will error if g and Wall have different sizes. 
"""
function Config(tf::Real, dt::Real, g::T; k=4.0e2::Real, walls=[]::Vector{Wall}) where {T}

	D::Int = length(g)

	TF = convert(Float64,tf)
	DT = convert(Float64,dt)
	K = convert(Float64,k)
	G = convert(SVector{D,Float64}, g)
	W = convert(Vector{Wall{D}}, unique(walls))


	Config(TF,DT,K,G,W)
end