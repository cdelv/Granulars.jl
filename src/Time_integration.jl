"""
This function does all the work. It's the core of the simulation.
- data: Array of Particles. 
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- vis_steps: How often to save a frame.
- file: Where to save the simulation information.
- save: Whether or not to save simulation data. 
"""
function Propagate(data::Vector{Particle}, conf::Config; vis_steps=2000::Int, file="Paraview/data"::String, save=false::Bool)
    
    # Time and Printing Variables
    t::Float64 = 0.0
    Print::Int = 0
    
    # Create Structure of Arrays
    particles = StructArray(data)
    
    # Cutoff For The Cell Lists.
    Cutoff::Float64 = 4*maximum(particles.rad)
    
    # Create Neighbor and Cell List
    system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
    list = neighborlist!(system) # Type Warning

    # Create Simetric Sparse Matrix for kundall spring forces calculation.
    # Sparse because most of the elements are zero and for large number of
    # particles we run out of memory.
    kundall_particles = ExtendableSparseMatrix(zeros(length(data), length(data)))
    kundall_walls = ExtendableSparseMatrix(zeros(length(data), length(conf.walls)))
    
    # Save Initial Condition. Save_step is Defined in Writte.jl
    if save
        Save_step(particles,file,Print,t); Print+=1
    end
    
    # Time Integration
    for i in 1:trunc(Int, conf.tf/conf.dt) # Number of Steps
        t+=conf.dt
        
        time_step(particles,conf,list,kundall_particles,kundall_walls)
        
        # Update Cell List
        update!(system, particles.r)
        list = neighborlist!(system)

        # Remove 0 entries from the sparse matrix.
        dropzeros!(kundall_particles)
        dropzeros!(kundall_walls)
        
        # Save Data Every vis_steps
        if i%vis_steps==0 && save 
            Save_step(particles,file,Print)
            Print+=1
        end
    end
    # Return The Final State Of The Simulation. 
    particles
end

"""
Performs one-time steps according to the PEFRL algorithm.
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- neighbor list: Neighbor list for a particle-to-particle interaction force calculation. 
- kundall_particles: Sparse symmetric matrix that stores the kundall spring distance for particle-particle interactions.
- kundall_walls: Sparse symmetric matrix that stores the kundall spring for particle-wall interactions.

The rotations algorithm is PEFRL:

- Optimized Forest–Ruth- and Suzuki-like algorithms for integration of motion in many-body systems, I.P. Omelyan, I.M. Mryglodab and R. Folk, 2002

PEFRL is O(4) in position and O(3) in velocity.

The rotations algorithm is:

- Algorithm for numerical integration of the rigid-body equations of motion, Igor P. Omelyan, 1998

This algortihm is O(3) in orientation and angular velocity and uses quaternions.
"""
function time_step(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    kundall_particles::ExtendableSparseMatrix{Float64, Int64},
    kundall_walls::ExtendableSparseMatrix{Float64, Int64})
    
    #PEFRL constants
    const1::Float64 = 0.1644986515575760     #ζ
    const3::Float64 = 0.1235692651138917e1   #χ
    const4::Float64 = -0.2094333910398989e-1 #λ
    const2::Float64 = (1-2*const4)/2         #(1-2λ)/2
    const5::Float64 = 1-2*(const3+const1);   #1-2*(ζ+χ)

    # Update Position. Move_r, Move_q, Move_w are defined in Particle.jl.
    for i in eachindex(particles)
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt,const1)
        @inbounds particles.w[i] = Move_w(particles.w[i],particles.τ[i],conf.dt,0.5)
        @inbounds particles.q[i] = Move_q(particles.q[i],particles.w[i],conf.dt)
    end
    # Calculate_Forces is defined in Forces.jl
    Calculate_Forces(particles,neighborlist,conf,kundall_particles,kundall_walls)
    # Update velocity. Move_v is defined in Particle.jl. 
    for i in eachindex(particles)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt,const2)
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt,const3)
    end
    Calculate_Forces(particles,neighborlist,conf,kundall_particles,kundall_walls)
    for i in eachindex(particles)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt,const4)
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt,const5)
    end
    Calculate_Forces(particles,neighborlist,conf,kundall_particles,kundall_walls)
    for i in eachindex(particles)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt,const4)
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt,const3)
    end
    Calculate_Forces(particles,neighborlist,conf,kundall_particles,kundall_walls)
    for i in eachindex(particles)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt,const2)
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt,const1)
        @inbounds particles.w[i] = Move_w(particles.w[i],particles.τ[i],conf.dt,0.5)
    end

    return nothing
end