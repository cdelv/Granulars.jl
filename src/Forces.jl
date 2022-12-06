"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
"""
function Calculate_Forces(particles::StructVector{Particle{D}}, neighborlist::Vector{Tuple{Int64, Int64, Float64}}, conf::Config{D}) where {D}
    
    @inbounds for i in eachindex(particles)
        # Reset Forces and Add Gravity.
        @inbounds particles.a[i] = conf.g

        # Calculate and Add Forces With Walls
        @inbounds particles.a[i] += Force_With_Walls(particles[i],conf)
    end

    # Calculate Forces Between Particles using the neighborlist.
    Force_With_Pairs(particles, conf, neighborlist)

    return nothing
end

"""
Calculate the force between pairs of particles using the neighbor list. 
- TO DO: IMPLEMENT FRICTION.
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
- i: an element of the neighbor list. 
"""
function Force_With_Pairs(particles::StructVector{Particle{D}}, conf::Config{D}, list::Vector{Tuple{Int64, Int64, Float64}}) where {D}
    
    for pair in list
        @inbounds i,j,d = pair
        @inbounds s = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that neighborlist hass a bigger cuttof. 
        if s > 0
            # Force to add.
            @inbounds F = Hertz_Force(s,conf)*normalize(particles.r[i]-particles.r[j])
            @inbounds F += Damping_Force(s,particles[i],particles[j],conf)

            #Newton 2 law. 
            @inbounds particles.a[i] += F/particles.m[i]
            @inbounds particles.a[j] -= F/particles.m[j]
        end
    end

    return nothing
end

"""
Uses the distance between a plane and a point to check for contact with the walls.
Then, it produces Hertz's force in the direction of the normal vector.
- TO DO: IMPLEMENT FRICTION AND DAMPING.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl.  
"""
function Force_With_Walls(particle::Particle{D}, conf::Config{D})::SVector{D, Float64} where {D}

	F::SVector{D, Float64} = zeros(SVector{D})

	for wall in conf.walls
		s::Float64 = particle.rad - norm(dot(particle.r-wall.Q,wall.n))
		if s > 0
			F += Hertz_Force(s,conf)*wall.n
            F += Damping_Force(s,particle,conf)
		end
	end
	F
end

"""
Hertz Force
- s: Interpenetration distance.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
"""
function Hertz_Force(s::Float64, conf::Config{D})::Float64 where {D}
	conf.K*s^1.5
end

"""
Damping Force in the collision of 2 particles. 
- s: Interpenetration distance.
- p1 and p2: particles interacting.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
function Damping_Force(s::Float64, p1::Particle{D}, p2::Particle{D}, conf::Config{D})::SVector{D, Float64} where {D}
    
    #Reduced Mass
    m12 = p1.m*p2.m/(p1.m+p2.m)

    #Relative Velocity
    V12 = p1.v-p2.v # carefull with angular velocity !

    -conf.gamma*sqrt(s)*m12*V12
end

"""
Damping Force in the collision of a particles and a wall. 
- s: Interpenetration distance.
- p1 : particles interacting with the wall.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
function Damping_Force(s::Float64, p1::Particle{D}, conf::Config{D})::SVector{D, Float64} where {D}
    -conf.gamma*sqrt(s)*p1.v # carefull with angular velocity !
end