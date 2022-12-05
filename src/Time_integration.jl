"""
This function does all the work. It's the core of the simulation.
- data: Array of Particles. 
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
- vis_steps: How often to save a frame.
- file: Where to save the simulation information.
- save: Whether or not to save simulation data. 
"""
function Propagate(data::Vector{Particle{D}}, conf::Config{D}; vis_steps=2000::Int, file="Paraview/data"::String, save=false::Bool) where {D}
    
    # Define Variables
    t::Float64 = 0.0 # Check if needed. Currently does nothing
    Print::Int = 0
    
    # Create Structure of Array
    particles = StructArray(data)
    
    # Calculate the Cutoff For The Cell Lists. MAKE A CALCULATION THAT MAKES SENSE
    Cutoff::Float64 = 6*maximum(particles.rad)
    
    # Create Neighbor and Cell List
    system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
    list = neighborlist!(system) # Type Warning
    
    # Save Initial Condition. Save_step is Deffined in Writte.jl
    if save
        Save_step(particles,file,Print); Print+=1
    end
    
    # Time Integration
    for i in 1:trunc(Int, conf.tf/conf.dt) # Number of Steps
        t+=conf.dt
        
        # Time Step
        PEFRL_time_step(particles,conf,list)
        
        # Update Cell List (Dont do it all steps)
        update!(system, particles.r)
        list = neighborlist!(system)
        
        # Save Data Every vis_steps
        if i%vis_steps==0 && save 
            Save_step(particles,file,Print); Print+=1
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
"""
function PEFRL_time_step(particles::StructVector{<:AbstractParticle{D}}, conf::Config{D}, neighborlist::Vector{Tuple{Int64, Int64, Float64}}) where {D}
    
    #PEFRL constants
    const1::Float64 = 0.1644986515575760     #ζ
    const3::Float64 = 0.1235692651138917e1   #χ
    const4::Float64 = -0.2094333910398989e-1 #λ
    const2::Float64 = (1-2*const4)/2         #(1-2λ)/2
    const5::Float64 = 1-2*(const3+const1);   #1-2*(ζ+χ)

    # Update Position. The map is faster than a For Loop for big systems. 
    # Also, it makes fewer allocations. Move_r is defined in Particle.jl 
    # and returns a SVector{D, Float64}.
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const1), LazyRows(particles))
    # Calculate_Forces is defined in Forces.jl
    Calculate_Forces(particles,neighborlist,conf)
    # Update velocity. Move_v is defined in Particle.jl. 
    # Also returns a SVector{D, Float64}.
    map(x -> x.v = Move_v(x.v,x.a,conf.dt,const2), LazyRows(particles))
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const3), LazyRows(particles))
    Calculate_Forces(particles,neighborlist,conf)
    map(x -> x.v = Move_v(x.v,x.a,conf.dt,const4), LazyRows(particles))
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const5), LazyRows(particles))
    Calculate_Forces(particles,neighborlist,conf)
    map(x -> x.v = Move_v(x.v,x.a,conf.dt,const4), LazyRows(particles))
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const3), LazyRows(particles))
	Calculate_Forces(particles,neighborlist,conf)
    map(x -> x.v = Move_v(x.v,x.a,conf.dt,const2), LazyRows(particles))
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const1), LazyRows(particles))

    return nothing
end