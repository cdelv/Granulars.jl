"""
This function does all the work. It's the core of the simulation.
- data: Array of Particles. 
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- vis_steps: How often to save a frame.
- file: Where to save the simulation information.
- save: Whether or not to save simulation data. 
- rot_seq: Rotation sequence to use for the orientation angle saving convention.
- beam_forces: Whether or not to create beams between overlaping particles at the begining. 
- fixed_spheres: List of indices of the spheres that wont move. They can rotate. 
- static: Make fixed spheres unable to rotate. 
"""
function Propagate(data::Vector{Particle}, 
    conf::Config; 
    vis_steps::Int=2000, 
    file::String="Paraview/data", 
    save::Bool=false, 
    rot_seq::Symbol=:XYZ, 
    beam_forces::Bool=false,
    fixed_spheres::Vector{Int64}=Int64[],
    static::Bool=false)::StructVector{<:AbstractParticle}
    
    # Time and Printing Variables.
    t::Float64 = 0.0
    Print::Int = 0
    
    # Create Structure of Arrays.
    particles = StructArray(data)
    
    # Cutoff For The Cell Lists.
    Cutoff::Float64 = 3.0*maximum(particles.rad)
    
    # Create Neighbor and Cell List.
    system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
    list = neighborlist!(system) # Type Warning!

    # Create dictionary for spring force calculation.
    friction_spring_particles = Dict{Tuple{Int64, Int64}, SVector{3, Float64}}()
    friction_spring_walls = Dict{Tuple{Int64, Int64}, SVector{3, Float64}}()

    # Dictionary that stores wich particles pairs have beam bonds.
    # Stores the index of the array of beams that correspond to the bond.
    beam_bonds = Dict{Tuple{Int64, Int64}, Int64}()

    # Stores the simulation beams
    beams = Beam[]

    # Creates beams between all ovelaping particles, Defined in Beams.jl
    if beam_forces
        Create_beams!(particles, list, beam_bonds, beams, conf)
    end

    # Stores the simulation beams (for performance)
    beams = StructArray(beams)
    
    # Save Initial Condition. Save_step is Defined in Writte.jl
    if save
        Save_step(particles,file,Print,t,rot_seq)
        Print+=1
    end

    # Needed for the LeapFrog Algorithm
    time_step_start!(particles,conf,list,friction_spring_particles,friction_spring_walls,beam_bonds,beams,fixed_spheres,static,t)
    
    # Time Integration.
    for i in 1:trunc(Int, conf.tf/conf.dt) # Number of Steps.
        t+=conf.dt
        
        Actions_Before_Time_Step!(particles,conf,list,friction_spring_particles,friction_spring_walls,beam_bonds,beams,fixed_spheres,static,t)

        time_step!(particles,conf,list,friction_spring_particles,friction_spring_walls,beam_bonds,beams,fixed_spheres,static,t)

        Actions_After_Time_Step!(particles,conf,list,friction_spring_particles,friction_spring_walls,beam_bonds,beams,fixed_spheres,static,t)

        # Update Cell List.
        update!(system, particles.r)
        list = neighborlist!(system)
        
        # Save Data Every vis_steps.
        if i%vis_steps==0 && save 
            Save_step(particles,file,Print,t,rot_seq) # Defined in Writte.jl
            Print+=1
        end
    end

    # Return The Final State Of The Simulation. 
    return particles
end

"""
Performs time-step initialization according to the LeapFrog algorithm.
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- neighborlist: Neighbor list for a particle-to-particle interaction force calculation. 
- friction_spring_particles: Dictionary that stores the spring distance for particle-particle interactions.
- friction_spring_walls: Dictionary that stores the spring distance for particle-wall interactions.

The traslation algorithm is LeapFrog, which is O(3) in position and velocity. 

The rotations algorithm is:

- Algorithm for numerical integration of the rigid-body equations of motion, Igor P. Omelyan, 1998

This algortihm is O(3) in orientation and angular velocity and uses quaternions. Similar to LeapFrog.
"""
function time_step_start!(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    friction_spring_particles::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    friction_spring_walls::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    beam_bonds::Dict{Tuple{Int64, Int64}, Int64},
    beams::StructVector{Beam},
    fixed_spheres::Vector{Int64},
    static::Bool,
    t::Float64)
    
    # Calculate_Forces is defined in Forces.jl
    Calculate_Forces!(particles,neighborlist,conf,friction_spring_particles,friction_spring_walls,beam_bonds,beams,t)

    # Remove forces and torques acting over static or fixed spheres
    for i in fixed_spheres
        particles.a[i] *= 0.0
        if static 
            particles.τ[i] *= 0.0
        end
    end

    # Initialize Velocity. Move_w, and Move_v are defined in Particle.jl.
    for i in eachindex(particles)
        @inbounds particles.w[i] = Move_w(particles.w[i],particles.τ[i],particles.I[i],conf.dt,-0.5)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt,-0.5)
    end

end


"""
Performs one time-step according to the LeapFrog algorithm.
- particles: StructArray of particles.
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- neighbor list: Neighbor list for a particle-to-particle interaction force calculation. 
- friction_spring_particles: Dictionary that stores the spring distance for particle-particle interactions.
- friction_spring_walls: Dictionary that stores the spring distance for particle-wall interactions.

The traslation algorithm is LeapFrog, which is O(3) in position and velocity. 

The rotations algorithm is:

- Algorithm for numerical integration of the rigid-body equations of motion, Igor P. Omelyan, 1998

This algortihm is O(3) in orientation and angular velocity and uses quaternions. Similar to LeapFrog.
"""
function time_step!(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    friction_spring_particles::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    friction_spring_walls::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    beam_bonds::Dict{Tuple{Int64, Int64}, Int64},
    beams::StructVector{Beam},
    fixed_spheres::Vector{Int64},
    static::Bool,
    t::Float64)

    #Calculate_Forces is defined in Forces.jl
    Calculate_Forces!(particles,neighborlist,conf,friction_spring_particles,friction_spring_walls,beam_bonds,beams,t)

    # Remove forces and torques acting over static or fixed spheres
    for i in fixed_spheres
        particles.a[i] *= 0.0
        if static 
            particles.τ[i] *= 0.0
        end
    end

    # Update Position and Velocity. Move_r, Move_q, Move_w, and Move_v are defined in Particle.jl.
    for i in eachindex(particles)
        @inbounds particles.w[i] = Move_w(particles.w[i],particles.τ[i],particles.I[i],conf.dt)
        @inbounds particles.v[i] = Move_v(particles.v[i],particles.a[i],conf.dt)
        
        @inbounds particles.r[i] = Move_r(particles.r[i],particles.v[i],conf.dt)
        @inbounds particles.q[i] = Move_q(particles.q[i],particles.w[i],conf.dt)
    end

    return nothing
end