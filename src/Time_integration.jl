"""
This function does all the work. Its the core of the simulation.
- data: Array of Particles. 
- conf: Simulation configuration, its a Conf struct. 
- vis_steps: How often to save a frame.
- file: Where to save the simulation information.
- save: Whether or not to save simulation data.  
"""
function Propagate(data::Vector{Particle{D}}, conf::Config{D}; vis_steps=2000::Int, file="Paraview/data"::String, save=false::Bool) where {D}
    
    # Define Variables
    t::Float64 = 0.0 # Check if needed, currently does nothing. 
    Print::Int = 0
    
    # Create Structure of Array
    particles = StructArray(data)
    
    # Calculate Cutoff for the Cell Lists. MAKE A CALCULATION THAT MAKES SENSE.
    Cutoff::Float64 = 6*maximum(particles.rad)
    
    # Create Neighbor and Cell List
    system = InPlaceNeighborList(x=particles.r, cutoff=Cutoff, parallel=false) # Explore Parallel Options
    list = neighborlist!(system) # Type Warning
    
    # Save Initial Condition. Save_step is deffined in Writte.jl
    if save
        Save_step(particles,file,Print); Print+=1
    end
    
    # Time Integration
    for i in 1:trunc(Int, conf.tf/conf.dt) # Number of steps
        t+=conf.dt
        
        # Time step. 
        PEFRL_time_step(particles,conf,list)
        
        # Update Cell List (Dont do it all steps)
        update!(system, particles.r)
        list = neighborlist!(system)
        
        # Save Data Every vis_steps
        if i%vis_steps==0 && save 
            Save_step(particles,file,Print); Print+=1
        end
    end
    # Return the final state of the simulation. 
    particles
end

"""
Performs 1 time steps according to the PEFRL algorithm.
- particles: StructArray of particles.
- conf: Simulation configuration, its a Conf struct. 
- neighborlist: Neighbor list for particle to particle interaction. 
"""
function PEFRL_time_step(particles::StructVector{<:AbstractParticle{D}}, conf::Config{D}, neighborlist::Vector{Tuple{Int64, Int64, Float64}}) where {D}
    #PEFRL constants
    const1::Float64 = 0.1644986515575760     #ζ
    const3::Float64 = 0.1235692651138917e1   #χ
    const4::Float64 = -0.2094333910398989e-1 #λ
    const2::Float64 = (1-2*const4)/2         #(1-2λ)/2
    const5::Float64 = 1-2*(const3+const1);   #1-2*(ζ+χ)

    # Update possition, map is faster than a for loop for big systems. 
    # Also it makes fewer allocations. Move_r is deffined in Particle.jl
    map(x -> x.r = Move_r(x.r,x.v,conf.dt,const1), LazyRows(particles))
    #Calculate_Forces is deffined in Forces.jl
    Calculate_Forces(particles,neighborlist,conf)
    # Update velocity. Move_v is deffined in Particle.jl
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