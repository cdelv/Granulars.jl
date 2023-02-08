"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- cundall_particles: Sparse symmetric matrix that stores the Cundall spring distance for particle-particle interactions.
- cundall_walls: Sparse symmetric matrix that stores the Cundall spring for particle-wall interactions.
- beam_bonds:
- beams:

COMPLETE
"""
function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, conf::Config,
    cundall_particles::ExtendableSparseMatrix{Float64, Int64},
    cundall_walls::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam})
    
    for i in eachindex(particles)
        # Calculate and Add Forces With Walls, torque and force is reset in the function
        Force_With_Walls(particles, i, conf, cundall_walls)
    end

    # Calculate Forces Between Particles using the neighborlist.
    Force_With_Pairs(particles, conf, neighborlist, cundall_particles,beam_bonds,beams)

    #particles.a[1] = 0*particles.a[1]
    #particles.τ[1] = 0*particles.τ[1]

    return nothing
end

"""
Calculate the force between pairs of particles using the neighbor list. 
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- cundall: Sparse symmetric matrix that stores the Cundall spring distance for particle-particle interactions.
- beam_bonds:
- beams:

COMPLETE
"""
function Force_With_Pairs(particles::StructVector{Particle}, conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    cundall::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam})
    
    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1],pair[2]) # For symetric acces to the Cundall distance matrix
        @inbounds j::Int64 = max(pair[1],pair[2]) # For symetric acces to the Cundall distance matrix
        @inbounds d::Float64 = pair[3]

        # Normal vector.
        @inbounds n::SVector{3, Float64} = unitary( particles.r[i] - particles.r[j] ) # The normal goes from j to i.

        # Calculate beam forces.
        if beam_bonds[i,j]!=0
            Beam_Force(particles, beams, i, j, beam_bonds[i,j], -n, conf) # -n because it has to go from i to j
            continue # No contact forces between beam bonded particles
        end

        # Interpenetration distance.
        # rᵢ - rⱼ means that goes from j to i. j=1 and i=2.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s < 0.0
            # Reset Cundall spring distance if there is no contact. 
            @inbounds cundall[i,j] = 0.0
            dropzeros!(cundall) # Remove 0 entries from the sparse matrix.
            continue
        end

        # Relative velocity. Carefull with angular velocity!
        @inbounds Vij::SVector{3, Float64} = (particles.v[i] + cross( Body_to_lab(particles.w[i],particles.q[i]), -particles.rad[i]*n ) 
            - (particles.v[j] + cross( Body_to_lab(particles.w[j],particles.q[j]), particles.rad[j]*n )))

        # Tangential vector and Reduced mass.
        t::SVector{3, Float64} = unitary( Vij - dot(Vij,n)*n )
        @inbounds mij::Float64 = particles.m[i]*particles.m[j]/(particles.m[i] + particles.m[j])

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,n)
        Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,mij,Vn,conf)
        if Fn < 0
            Fn = 0.0
        end

        # Calculate tangencial forces (Cundall friction force)
        Vt::Float64 = dot(Vij,t)
        @inbounds cundall[i,j] += Vt*conf.dt # Add distance to Cundall spring
        @inbounds Ft::Float64 = Cundall_friction(cundall[i,j], Fn, conf)

        # Total force
        F::SVector{3, Float64} = Fn*n + Ft*t

        # Add force (Newton 2 law)
        @inbounds particles.a[i] += F/particles.m[i]
        @inbounds particles.a[j] -= F/particles.m[j]

        # Add torque
        @inbounds particles.τ[i] += Lab_to_body(cross(-particles.rad[i]*n, F), particles.q[i])
        @inbounds particles.τ[j] += Lab_to_body(cross( particles.rad[j]*n,-F), particles.q[j])
    end

    return nothing
end

"""
Uses the distance between a plane and a point to check for contact with the walls.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
- cundall: Sparse symmetric matrix that stores the Cundall spring for particle-wall interactions. 
"""
function Force_With_Walls(particles::StructVector{Particle}, i::Int64, conf::Config,
    cundall::ExtendableSparseMatrix{Float64, Int64})
    
    # Reset torques and set gravity to reset forces.
    @inbounds F::SVector{3, Float64} = conf.g*particles.m[i]
    T::SVector{3, Float64} = zeros(SVector{3})

    for j in eachindex(conf.walls)

        # Interpenetration distance.
        @inbounds s::Float64 = particles.rad[i] - dot(particles.r[i]-conf.walls[j].Q, conf.walls[j].n)

        # Check for contact.
        if s < 0.0
            # Reset Cundall spring distance if theres no contact. 
            cundall[i,j] = 0.0
            dropzeros!(cundall) # Remove 0 entries from the sparse matrix.
            continue
        end

        # Relative velocity. Carefull with angular velocity! The minus sing is due to the direction of the normal.
        @inbounds Vij::SVector{3, Float64} = particles.v[i] + cross(Body_to_lab(particles.w[i],particles.q[i]), -particles.rad[i]*conf.walls[j].n)

        # Tangential vector. The normal one is conf.walls[j].n, it enters the particle.
        t::SVector{3, Float64} = unitary(Vij - dot(Vij,conf.walls[j].n)*conf.walls[j].n)

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,conf.walls[j].n)
        @inbounds Fn::Float64 = Hertz_Force(s,conf) + Damping_Force(s,particles.m[i],Vn,conf) # Reduced mass is m (wall with infinite mass).
        if Fn < 0
            Fn = 0.0
        end

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        @inbounds cundall[i,j] += Vt*conf.dt # Add distance to Cundall spring
        @inbounds Ft::Float64 = Cundall_friction(cundall[i,j], Fn, conf)
        
        # Add Total force and Torque of this wall
        @inbounds F += Fn*conf.walls[j].n + Ft*t
        @inbounds T += cross(-particles.rad[i]*conf.walls[j].n, Fn*conf.walls[j].n+Ft*t)
    end

    @inbounds particles.a[i] = F/particles.m[i]
    @inbounds particles.τ[i] = Lab_to_body(T,particles.q[i])
    
    return nothing
end

"""
Hertz Force.
- s: Interpenetration distance.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
"""
@inline function Hertz_Force(s::Float64, conf::Config)::Float64
	conf.K*s*sqrt(s) # s^1.5 but faster
end

"""
Damping Force. 
- s: Interpenetration distance.
- mij: Reduced mass.
- Vn: Normal component of relative velocity. 
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
@inline function Damping_Force(s::Float64, mij::Float64, Vn::Float64, conf::Config)::Float64
    -conf.gamma*sqrt(s)*mij*Vn
end

"""
Calculates the friction force, wheder it is cinetic or static using the Cundall friction model. 
- X_cundall: Cundall spring elongation. 
- Fn: Normal force aplied to the particle.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
@inline function Cundall_friction(X_cundall::Float64, Fn::Float64, conf::Config)::Float64
    Ft::Float64 = -conf.K_cundall*X_cundall
    Ftmax::Float64 = conf.mu*Fn         

    # Check if friction is static or dynamic. 
    # TO DO: Diferenciate between the 2 coeficients BETTER!
    if abs(Ft) > Ftmax 
        Ft = sign(Ft)*Ftmax
    end
    Ft
end