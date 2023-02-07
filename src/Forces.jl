"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- cundall_particles: Sparse symmetric matrix that stores the Cundall spring distance for particle-particle interactions.
- cundall_walls: Sparse symmetric matrix that stores the Cundall spring for particle-wall interactions.
"""
function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    conf::Config,
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

    return nothing
end

"""
Calculate the force between pairs of particles using the neighbor list. 
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- cundall: Sparse symmetric matrix that stores the Cundall spring distance for particle-particle interactions.
"""
function Force_With_Pairs(particles::StructVector{Particle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    cundall::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam})
    
    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1],pair[2]) # For symetric acces to the Cundall distance matrix
        @inbounds j::Int64 = max(pair[1],pair[2]) # For symetric acces to the Cundall distance matrix
        @inbounds d::Float64 = pair[3]

        # Interpenetration distance.
        # rᵢ - rⱼ means that goes from j to i. j=1 and i=2.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Normal vector.
        @inbounds n::SVector{3, Float64} = unitary( particles.r[i] - particles.r[j] ) # The normal goes from j to i.

        # Calculate beam forces.
        if beam_bonds[i,j]!=0
            Beam_Force(particles, beams, i, j, beam_bonds[i,j], -n, conf) # -n because it has to go from i to j
            continue # No contact forces between beam bonded particles
        end

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
        @inbounds particles.τ[i] += Lab_to_body(cross(-particles.rad[i]*n, F), particles.q[i]) #Check this SIGNS!!
        @inbounds particles.τ[j] += Lab_to_body(cross(particles.rad[j]*n, -F), particles.q[j])
    end

    particles.τ[1] = SVector(0.0,0.0,0.0)
    particles.a[1] = SVector(0.0,0.0,0.0)

    particles.τ[100] = SVector(0.0,0.0,0.0)
    particles.a[100] = SVector(0.0,0.0,0.0)

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
        @inbounds F += Fn*conf.walls[j].n

        # Calculate tangencial forces (Kundal friction force)
        Vt::Float64 = dot(Vij,t)
        cundall[i,j] += Vt*conf.dt # Add distance to Cundall spring
        @inbounds Ft::Float64 = Cundall_friction(cundall[i,j], Fn, conf)
        F += Ft*t
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


"""
- particles: StructArray of particles.
- i:
- j:
- k:
- n:
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
"""
function Beam_Force(particles::StructVector{Particle}, 
    beams::StructVector{Beam}, i::Int64, j::Int64, k::Int64, 
    n::SVector{3,Float64}, conf::Config)

    # Calculate displacements in the beam frame.
    dxi::SVector{3,Float64} = Lab_to_Beam(n, beams.r_i[k] - particles.r[i])
    dΩi::SVector{3,Float64} = Lab_to_Beam(n, quat_to_angle(beams.q_i[k]) - quat_to_angle(particles.q[i]))
    dxj::SVector{3,Float64} = Lab_to_Beam(n, beams.r_j[k] - particles.r[j])
    dΩj::SVector{3,Float64} = Lab_to_Beam(n, quat_to_angle(beams.q_j[k]) - quat_to_angle(particles.q[j]))

    # Create the 12x12 transformation matrix and calculate forces and torques.
    @inbounds Δs::SVector{12,Float64} = SVector(dxi[1], dxi[2], dxi[3], dΩi[1], dΩi[2], dΩi[3], dxj[1], dxj[2], dxj[3], dΩj[1], dΩj[2], dΩj[3])
    F::SVector{12,Float64} = K_beam(beams.A[k], beams.L[k], conf.E, conf.G)*Δs #cte k matrix, maybe compute once

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Beam_to_Lab(n, SVector(F[1], F[2], F[3]))/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Beam_to_Lab(n, SVector(F[4], F[5], F[6])), particles.q[i])

    @inbounds particles.a[j] += Beam_to_Lab(n, SVector(F[7], F[8], F[9]))/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Beam_to_Lab(n, SVector(F[10],F[11],F[12])), particles.q[j])

    γ = 0.1
    @inbounds particles.a[i] += -γ*particles.v[i]
    @inbounds particles.τ[i] += -γ*particles.w[i]

    @inbounds particles.a[j] += -γ*particles.v[j]
    @inbounds particles.τ[j] += -γ*particles.w[j]

    # Update the beam information (NOT NEEDED DUE TO THE WAY THE K MATRIX WORKS)
    #beams.r_i[k] = particles.r[i]
    #beams.r_j[k] = particles.r[j]
    #beams.q_i[k] = particles.q[i]
    #beams.q_j[k] = particles.q[j]

    return nothing
end