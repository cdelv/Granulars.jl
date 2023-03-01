include("../../src/Granulars.jl")

function main()
    # Simulation parameters
    g = [0.0,-10.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx,Ly,Lz, E=1e9, ν=-0.5)

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0, E=1e7)
    p2 = Particle(r=[13, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0, E=1e6)
    p3 = Particle(r=[11.5, 13, 10], v=[10.0,0.0,0], w=[0,0,0], rad=2.0, E=1e6)
    particles = [p1, p2, p3]

    # Estimate a good time step
    vis_steps = 900
    frames = 200
    dt = 0.7*PWaveTimeStep(particles)
    println("dt = ", dt)
    
    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, walls=walls, mu=0.4, beam_forces=true, thorsten_damping=true, en=0.9, v=0.1)

    # Run the simulation
    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();