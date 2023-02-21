"""
Main force calculation routine. It calls all the different force interactions. 
- particles: StructArray of particles.
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- mindlin_F_particles: Sparse symmetric matrix that stores the Mindlin spring distance for particle-particle interactions.
- mindlin_F_walls: Sparse symmetric matrix that stores the Mindlin spring for particle-wall interactions.
- beam_bonds: Sparse symmetric matrix that stores wich beam connects with each particle pair.
- beams: StructArray of beams between particles.
"""
function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, conf::Config,
    mindlin_F_particles::ExtendableSparseMatrix{Float64, Int64},
    mindlin_F_walls::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam},
    t::Float64)
    
    # Reset force acting upon the wall
    for i in eachindex(conf.walls)
        conf.walls[i] = Set_F(conf.walls[i], zeros(SVector{3}))
    end
    
    for i in eachindex(particles)
        # Calculate and Add Forces With Walls, torque and force is reset in the function
        Force_With_Walls(particles, i, conf, mindlin_F_walls)
    end

    # Calculate Forces Between Particles using the neighborlist.
    Force_With_Pairs(particles, conf, neighborlist, mindlin_F_particles,beam_bonds,beams)

    return nothing
end

"""
Calculate the force between pairs of particles using the neighbor list. 
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.  
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- mindlin_F: Sparse symmetric matrix that stores the Mindlin spring distance for particle-particle interactions.
- beam_bonds: Sparse symmetric matrix that stores wich beam connects with each particle pair.
- beams: StructArray of beams between particles.
"""
function Force_With_Pairs(particles::StructVector{Particle}, conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    mindlin_F::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam})
    
    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1],pair[2]) # For symetric acces to the Mindlin Force matrix
        @inbounds j::Int64 = max(pair[1],pair[2]) # For symetric acces to the Mindlin Force matrix
        @inbounds d::Float64 = pair[3]

        # Calculate beam forces.
        if beam_bonds[i,j]!=0
            Beam_Force(particles, beams, i, j, beam_bonds[i,j]) #Defined in Beams.jl
            continue # No contact forces between beam bonded particles
        end

        # Interpenetration distance.
        # rᵢ - rⱼ means that goes from j to i. j=1 and i=2.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s <= 0.0
            # Reset Cundall spring distance if there is no contact. 
            @inbounds mindlin_F[i,j] = 0.0
            dropzeros!(mindlin_F) # Remove 0 entries from the sparse matrix.
            continue
        end

        # Antiracheting 
        # https://gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/ScGeom.cpp#L65-94
        # proposed by McNamara and co-workers. dicusion -> DOI 10.1103/PhysRevE.77.031304
        α::Float64 = (particles.rad[i] + particles.rad[j])/(particles.rad[i] + particles.rad[j] - s)

        # Normal vector.
        @inbounds n::SVector{3, Float64} = unitary(particles.r[i] - particles.r[j]) # The normal goes from j to i.

        # Relative velocity. Carefull with angular velocity!
        @inbounds Vij::SVector{3, Float64} = (α*particles.v[i] + cross( Body_to_lab(particles.w[i],particles.q[i]), -particles.rad[i]*n ) 
                    - (α*particles.v[j] + cross( Body_to_lab(particles.w[j],particles.q[j]), particles.rad[j]*n )))

        # Reduced mass, radius, young modulus, and shear modulus
        @inbounds mij::Float64 = particles.m[i]*particles.m[j]/(particles.m[i] + particles.m[j])
        @inbounds Rij::Float64 = particles.rad[i]*particles.rad[j]/(particles.rad[i] + particles.rad[j])
        @inbounds Eij::Float64 = particles.E[i]*particles.E[j]/((1-particles.ν[i]^2)*particles.E[j]+(1-particles.ν[j]^2)*particles.E[i])
        @inbounds Gij::Float64 = particles.G[i]*particles.G[j]/((2-particles.ν[i])*particles.G[j]+(2-particles.ν[j])*particles.G[i])

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,n)
        Fn::Float64 = max(0.0, Hertz_Force(s,Eij,Rij) - Normal_Damping_Force(s, mij, Eij, Rij, Vn, conf.en))

        # Calculate tangencial forces
        Vt::SVector{3, Float64} = Vij - dot(Vij,n)*n
        Ft::SVector{3, Float64} = Mindlin_Shear_And_friction(mindlin_F, i, j, Vt, Gij, Rij, s, mij, Fn, conf)

        # Total force
        F::SVector{3, Float64} = Fn*n + Ft

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
- i: Index of the particle on wich the force is being calculated.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
- mindlin_F: Sparse symmetric matrix that stores the Mindlin spring for particle-wall interactions. 
"""
function Force_With_Walls(particles::StructVector{Particle}, i::Int64, conf::Config,
    mindlin_F::ExtendableSparseMatrix{Float64, Int64})
    
    # Reset torques and forces.
    @inbounds F::SVector{3, Float64} = zeros(SVector{3})
    T::SVector{3, Float64} = zeros(SVector{3})

    for j in eachindex(conf.walls)

        # Interpenetration distance.
        @inbounds s::Float64 = particles.rad[i] - dot(particles.r[i]-conf.walls[j].Q, conf.walls[j].n)

        # Check for contact.
        if s <= 0.0
            # Reset Cundall spring distance if theres no contact. 
            mindlin_F[i,j] = 0.0
            dropzeros!(mindlin_F) # Remove 0 entries from the sparse matrix.
            continue
        end

        # Relative velocity. Carefull with angular velocity! The minus sing is due to the direction of the normal.
        @inbounds Vij::SVector{3, Float64} = particles.v[i] - cross(Body_to_lab(particles.w[i],particles.q[i]), particles.rad[i]*conf.walls[j].n)

        # Reduced mass, radius, young modulus, and shear modulus
        @inbounds mij::Float64 = particles.m[i] # Reduced mass is m (wall with infinite mass).
        @inbounds Rij::Float64 = particles.rad[i] # Reduced radius is rad (wall with infinite radius).
        @inbounds Eij::Float64 = particles.E[i]*conf.walls[j].E/((1-particles.ν[i]^2)*conf.walls[j].E+(1-conf.walls[j].ν^2)*particles.E[i])
        @inbounds Gij::Float64 = particles.G[i]*conf.walls[j].G/((2-particles.ν[i])*conf.walls[j].G+(2-conf.walls[j].ν)*particles.G[i])

        # Calculate normal forces.
        Vn::Float64 = dot(Vij,conf.walls[j].n)
        Fn::Float64 = max(0.0, Hertz_Force(s,Eij,Rij) - Normal_Damping_Force(s, mij, Eij, Rij, Vn, conf.en))

        # Calculate tangencial forces
        Vt::SVector{3, Float64} = Vij - Vn*conf.walls[j].n
        Ft::SVector{3, Float64} = Mindlin_Shear_And_friction(mindlin_F, i, j, Vt, Gij, Rij, s, mij, Fn, conf)
        
        # Add Total force and Torque of this wall
        @inbounds F += Fn*conf.walls[j].n + Ft
        @inbounds T += cross(-particles.rad[i]*conf.walls[j].n, Fn*conf.walls[j].n + Ft)

        # Add force acting on the wall
        conf.walls[j] = Set_F(conf.walls[j], conf.walls[j].F - Fn*conf.walls[j].n - Ft)
    end

    # Update force acting on particle
    @inbounds particles.a[i] = conf.g + F/particles.m[i]
    @inbounds particles.τ[i] = Lab_to_body(T,particles.q[i])
    
    return nothing
end

"""
Hertz Force.
- s: Interpenetration distance.
- Eij: Reduced young modulus.
- Rij: Reduced radius.  
"""
@inline function Hertz_Force(s::Float64, Eij::Float64, Rij::Float64)::Float64
    4.0*Eij*sqrt(Rij*s)*s/3.0 # s^1.5 but faster
end

"""
Normal Damping Force. 

Taken from the book Granular Dynamics, Contact Mechanics and Particle System Simulations by Colin Thornton
equations (4.11 - 4.15)

- s: Interpenetration distance.
- mij: Reduced mass.
- Eij: Reduced young modulus.
- Rij: Reduced radius.
- Vn: Normal component of relative velocity. 
- en: Restitution coeficient.
"""
@inline function Normal_Damping_Force(s::Float64, mij::Float64, Eij::Float64, Rij::Float64, Vn::Float64, en::Float64)::Float64
    2.0*γ(en)*sqrt(mij*2.0*Eij*sqrt(Rij*s))*Vn
end

"""
Calculates the mindlin elastic shear force and kinetic friction.

Taken from the book Granular Dynamics, Contact Mechanics and Particle System Simulations by Colin Thornton
equations (4.11 - 4.15)

Not doing the checks for the normal force evolution...

- mindlin_F: Sparse simetric matrix for the Mindlin force acumulation.
- i: Index of the particle on wich the force is being calculated.
- j: Index of the other particle or wall on wich the force is being calculated.
- Vt: Tangential velocity.
- Gij: Reduced shear modulus.
- Rij: Reduced radius
- mij: Reduced mass.
- s: Interpenetration distance.
- Fn: Normal force magnitude. 
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl.

Shpuld need to rotate previus shear force but not duing it. See:
https://gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/ScGeom.cpp#L15-23
https://gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/ScGeom.cpp#L37-62
https://gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/HertzMindlin.cpp#L363-382
"""
@inline function Mindlin_Shear_And_friction(
    mindlin_F::ExtendableSparseMatrix{Float64, Int64},
    i::Int64, 
    j::Int64, 
    Vt::SVector{3, Float64}, 
    Gij::Float64, 
    Rij::Float64,
    mij::Float64, 
    s::Float64,
    Fn::Float64,
    conf::Config)::SVector{3, Float64}

    kt::Float64 = 8.0*Gij*sqrt(Rij*s)

    @inbounds mindlin_F[i,j] -= conf.dt
    @inbounds Ft::SVector{3, Float64} = (mindlin_F[i,j]*kt - 2.0*γ(conf.en)*sqrt(mij*kt))*Vt

    # TO DO: Diferenciate between the 2 friction coeficients!
    if norm(Ft) > Fn*conf.mu 
        return -Fn*conf.mu*Vt/norm(Vt) 
    end
    Ft
end