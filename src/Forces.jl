"""
Main force calculation ruttine. It calls oall the different force interactions. 
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct. 
- neighborlist: Neighbor list for particle to particle interaction. 
"""
function Calculate_Forces(particles::StructVector{Particle{D}}, neighborlist::Vector{Tuple{Int64, Int64, Float64}}, conf::Config{D}) where {D}
    # Reset Forces and Add Gravity. The map is much faster than the foor loop and avoids allocations. 
    map(x -> x.a = conf.g, LazyRows(particles)) 

    # Calculate and Add Forces With Walls
    map(x -> x.a += Force_With_Walls(x,conf), LazyRows(particles))

    # Calculate Forces Between Particles using the neighborlist. 
    map(x -> Force_With_Pairs(particles, conf, x), neighborlist)
end

"""
Calculate the force between pairs of particles using the neighborlist. 
TO DO: IMPLEMENT FRICTION AND DAMPING.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct. 
- i: an element of the neighborlist. 
"""
function Force_With_Pairs(particles::StructVector{Particle{D}}, conf::Config{D}, i::Tuple{Int64, Int64, Float64}) where {D}
    
    s = LazyRow(particles, i[1]).rad + LazyRow(particles, i[2]).rad - i[3]
    
    # Check for contact. Remember that neighborlist hass a bigger cuttof. 
    if s > 0
    	# Force to add. 
    	F = Hertz_Force(s,conf)*normalize(LazyRow(particles, i[1]).r-LazyRow(particles, i[2]).r)
        #Newton 2 law. 
        LazyRow(particles, i[1]).a += F
        LazyRow(particles, i[2]).a -= F
    end
    return nothing
end

"""
Uses the distance between a plain and a point to check for contact.
Then, it aplies a normal Hertz force to the contact.
TO DO: IMPLEMENT FRICTION AND DAMPING.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct. 
"""
function Force_With_Walls(particle::LazyRow{Particle{D}}, conf::Config{D})::SVector{D, Float64} where {D}

	F::SVector{D, Float64} = zeros(SVector{D})

	for wall in conf.walls
		s::Float64 = particle.rad - norm(dot(particle.r-wall.Q,wall.n))
		if s > 0
			F += Hertz_Force(s,conf)*wall.n
		end
	end
	F
end

"""
Hertz Force
- s: interpenetration distance
- conf: Simulation configuration, its a Conf struct. 
"""
function Hertz_Force(s::Float64, conf::Config{D}) where {D}
	conf.K*s^1.5
end