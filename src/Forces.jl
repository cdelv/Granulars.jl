"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- kundall_particles: Sparse symmetric matrix that stores the kundall spring distance for particle-particle interactions.
- kundall_walls: Sparse symmetric matrix that stores the kundall spring for particle-wall interactions.
"""
function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    conf::Config,
    kundall_particles::ExtendableSparseMatrix{Float64, Int64},
    kundall_walls::ExtendableSparseMatrix{Float64, Int64})
    
    for i in eachindex(particles)
        # Calculate and Add Forces With Walls, torque and force is reset in the function
        Force_With_Walls(particles, i, conf, kundall_walls)
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
function Force_With_Pairs(particles::StructVector{Particle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    kundall::ExtendableSparseMatrix{Float64, Int64})
    
    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1],pair[2]) # For symetric acces to the kundall distance matrix
        @inbounds j::Int64 = max(pair[1],pair[2]) # For symetric acces to the kundall distance matrix
        @inbounds d::Float64 = pair[3]

        # Interpenetration distance
        # rᵢ - rⱼ means that goes from j to i.  # j=1 and i=2 
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s < 0.0
            # Reset Kundall spring distance if theres no contact. 
            @inbounds kundall[i,j] = 0.0
            continue
        end

        # Normal and Relative velocity. Carefull with angular velocity!
        @inbounds n::SVector{3, Float64} = unitary( particles.r[i] - particles.r[j] ) # The minus sing is due to the direction of the normal j to i
        #                                                                   to account for the contact point (r-s/2)
        @inbounds Vij::SVector{3, Float64} = particles.v[i] + cross( particles.w[i], -(particles.rad[i]-s/2)*n ) - (particles.v[j] + cross( particles.w[j], (particles.rad[j]-s/2)*n ))

        # Tangential vector and Reduced mass
        t::SVector{3, Float64} = unitary( Vij - dot(Vij,n)*n )
        @inbounds mij::Float64 = particles.m[i]*particles.m[j]/(particles.m[i] + particles.m[j])

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,n)
        Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,mij,Vn,conf)

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        @inbounds kundall[i,j] += Vt*conf.dt*0.25 # Add distance to Kundall spring
        @inbounds Ft::Float64 = Kundall_friction(kundall[i,j], Fn, conf)

        # Total force
        F::SVector{3, Float64} = Fn*n + Ft*t

        # Add force (Newton 2 law) check sings
        @inbounds particles.a[i] += F/particles.m[i]
        @inbounds particles.a[j] -= F/particles.m[j]

        # Add torque
        @inbounds particles.τ[i] += cross(-(particles.rad[i]-s/2)*n, F)
        @inbounds particles.τ[j] += cross( (particles.rad[j]-s/2)*n, -F)
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
function Force_With_Walls(particles::StructVector{Particle}, i::Int64, conf::Config,
    kundall::ExtendableSparseMatrix{Float64, Int64})
    
    # Reset torques and set gravity to reset forces
    @inbounds F::SVector{3, Float64} = conf.g*particles.m[i]
    T::SVector{3, Float64} = zeros(SVector{3})

    for j in eachindex(conf.walls)

        # Interpenetration distance
        @inbounds s::Float64 = particles.rad[i] - dot(particles.r[i]-conf.walls[j].Q, conf.walls[j].n)

        # Check for contact
        if s < 0.0
            # Reset Kundall spring distance if theres no contact. 
            kundall[i,j] = 0.0
            continue
        end

        # Reduced mass and relative velocity. Carefull with angular velocity! The minus sing is due to the direction of the normal
        @inbounds Vij::SVector{3, Float64} = particles.v[i] + cross(particles.w[i], -(particles.rad[i]-s/2)*conf.walls[j].n)

        # Tangential vector. The normal one is conf.walls[j].n it enters the particle
        t::SVector{3, Float64} = unitary( Vij - dot(Vij,conf.walls[j].n)*conf.walls[j].n )

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,conf.walls[j].n)
        @inbounds Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,particles.m[i],Vn,conf) # Reduced mass is m (wall with infinite mass).
        @inbounds F += Fn*conf.walls[j].n
        if Fn < 0
            Fn = 0.0
        end

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        kundall[i,j] += Vt*conf.dt*0.25 # Add distance to Kundall spring
        @inbounds Ft::Float64 = Kundall_friction(kundall[i,j], abs(dot(F,conf.walls[j].n)), conf)
        F += Ft*t
        @inbounds T += cross(-(particles.rad[i]-s/2)*conf.walls[j].n, F)
    end

    @inbounds particles.a[i] = F/particles.m[i]
    @inbounds particles.τ[i] = T
    
    return nothing
end

"""
Hertz Force
- s: Interpenetration distance.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
"""
@inline function Hertz_Force(s::Float64, conf::Config)::Float64
	conf.K*s*sqrt(s) # ^1.5 but faster
end

"""
Damping Force in the collision of 2 particles. 
- s: Interpenetration distance.
- p1 and p2: particles interacting.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
@inline function Damping_Force(s::Float64, mij::Float64, Vn::Float64, conf::Config)::Float64
    -conf.gamma*sqrt(s)*mij*Vn
end

"""
Calculates the friction force, wheder it is cinetic or static using the Kundall force. 
- kundallX: Kundall spring elongation. 
- Fn: Normal force aplied to the particle.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
@inline function Kundall_friction(kundallX::Float64, Fn::Float64, conf::Config)::Float64
    Ft::Float64 = -conf.K_kundall*kundallX
    Ftmax::Float64 = conf.mu*abs(Fn)         

    # Check if friction is static or dynamic. TO DO: Diferenciate between the 2 coeficients.
    if abs(Ft) > Ftmax 
        Ft = sign(Ft)*Ftmax
    end
    Ft
end