include("../../src/Granulars.jl")

function main()
    # Simulation parameters
    g = [0.0,0.0,0.0]
    
    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,0.0,0], w=[1.0,0,0], rad=2.0, E=1e8, ν=0.2)
    p2 = Particle(r=[13, 10, 10], v=[0.0,0.0,0], w=[1.0,0,0], rad=2.0, E=1e8, ν=0.2)
    particles = [p1, p2]

    walls = Create_Box(25,25,25)

    # Estimate a good time step
    vis_steps = 1000
    frames = 200
    dt = 0.00001
    println("dt = ", dt)
    
    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, walls=walls, beam_forces=true, beam_damping=true)

    # Run the simulation
    Propagate([p1,p2], conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();