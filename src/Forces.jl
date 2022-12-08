"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- kundall_particles: Sparse symmetric matrix that stores the kundall spring distance for particle-particle interactions.
- kundall_walls: Sparse symmetric matrix that stores the kundall spring for particle-wall interactions.
"""
@inline function Calculate_Forces(particles::StructVector{Particle{D}}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    conf::Config{D},
    kundall_particles::ExtendableSparseMatrix{Float64, Int64},
    kundall_walls::ExtendableSparseMatrix{Float64, Int64}) where {D}
    
    for i in eachindex(particles)
        # Reset Forces and Add Gravity.
        particles.a[i] = conf.g

        # Calculate and Add Forces With Walls
        particles.a[i] += Force_With_Walls(particles[i], i, conf, kundall_walls)
    end

    # Calculate Forces Between Particles using the neighborlist.
    Force_With_Pairs(particles, conf, neighborlist, kundall_particles)

    return nothing
end

"""
Calculate the force between pairs of particles using the neighbor list. 
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- kundall: Sparse symmetric matrix that stores the kundall spring distance for particle-particle interactions.
"""
function Force_With_Pairs(particles::StructVector{Particle{D}}, 
    conf::Config{D}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    kundall::ExtendableSparseMatrix{Float64, Int64}) where {D}
    
    for pair in neighborlist
        i::Int64 = min(pair[1],pair[2]) # For symetric acces to the kundall distance matrix
        j::Int64 = max(pair[1],pair[2]) # For symetric acces to the kundall distance matrix
        d::Float64 = pair[3]

        # Interpenetration distance
        s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that neighborlist hass a bigger cuttof. 
        if s < 0.0
            # Reset Kundall spring distance if theres no contact. 
            kundall[i,j] = 0.0
            continue
        end

        # Reduced mass and relative velocity. Carefull with angular velocity!
        Vij::SVector{D, Float64} = particles.v[i] - particles.v[j]
        mij::Float64 = particles.m[i]*particles.m[j]/(particles.m[i] + particles.m[j])

        # Normal and tangential vectors
        n::SVector{D, Float64} = normalize( particles.r[i] - particles.r[j] )
        t::SVector{D, Float64} = normalize( Vij - dot(Vij,n)*n )

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,n)
        Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,mij,Vn,conf)
        if Fn < 0
            Fn = 0
        end

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        kundall[i,j] += Vt*conf.dt # Add distance to Kundall spring
        Ft = Kundall_friction(kundall[i,j], Fn, conf) + Damping_Force(s,mij,Vt,conf) # IS THIS CORRECT?

        # Total force
        F::SVector{D, Float64} = Fn*n + Ft*t

        # Add force (Newton 2 law) 
        particles.a[i] += F/particles.m[i]
        particles.a[j] -= F/particles.m[j]

        # Add torque
    end

    return nothing
end

"""
Uses the distance between a plane and a point to check for contact with the walls.
Then, it produces Hertz's force in the direction of the normal vector.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
- kundall: Sparse symmetric matrix that stores the kundall spring for particle-wall interactions. 
"""
function Force_With_Walls(particle::Particle{D}, i::Int64, conf::Config{D},
    kundall::ExtendableSparseMatrix{Float64, Int64})::SVector{D, Float64} where {D}

    F::SVector{D, Float64} = zeros(SVector{D})

    for j in eachindex(conf.walls)

        # Interpenetration distance
        s::Float64 = particle.rad - norm(dot(particle.r-conf.walls[j].Q,conf.walls[j].n))

        # Check for contact
        if s < 0.0
            # Reset Kundall spring distance if theres no contact. 
            kundall[i,j] = 0.0
            continue
        end

        # Reduced mass and relative velocity. Carefull with angular velocity!
        Vij::SVector{D, Float64} = particle.v

        # Tangential vector. The normal one is conf.walls[j].n
        t::SVector{D, Float64} = normalize( Vij - dot(Vij,conf.walls[j].n)*conf.walls[j].n )

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,conf.walls[j].n)
        Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,particle.m,Vn,conf)# Reduce mass is m (wall with infinite mass).
        if Fn < 0
            Fn = 0.0
        end

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        kundall[i,j] += Vt*conf.dt # Add distance to Kundall spring   # IS THIS CORRECT?
        Ft::Float64 = Kundall_friction(kundall[i,j], Fn, conf) + Damping_Force(s,particle.m,Vt,conf) 

        @inbounds F += Fn*conf.walls[j].n + Ft*t
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
function Damping_Force(s::Float64, mij::Float64, Vn::Float64, conf::Config{D})::Float64 where {D}
    -conf.gamma*sqrt(s)*mij*Vn
end

"""
Calculates the friction force, wheder it is cinetic or static using the Kundall force. 
- kundallX: Kundall spring elongation. 
- Fn: Normal force aplied to the particle.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
function Kundall_friction(kundallX::Float64, Fn::Float64, conf::Config{D})::Float64 where {D}
    Ft::Float64 = -conf.K_kundall*kundallX
    Ftmax::Float64 = conf.mu*abs(Fn)         

    # Check if friction is static or dynamic. TO DO: Diferenciate between the 2 coeficients.
    if abs(Ft) > Ftmax 
        Ft = sign(Ft)*Ftmax
    end
    Ft
end