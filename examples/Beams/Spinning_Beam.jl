include("../../src/Granulars.jl")

function main()
    # Simulation parameters
    g = [0.0,0.0,0.0]

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,-10.0,0], w=[0,0,0], rad=2.0, E=100000, ν=-0.5)
    p2 = Particle(r=[13, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0, E=100000, ν=-0.5)
    p3 = Particle(r=[16, 10, 10], v=[0.0,10.0,0], w=[0,0,0], rad=2.0, E=100000, ν=-0.5)
    particles = [p1,p2,p3]

    dt = 0.0001
    frames = 200
    vis_steps = 200

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, beam_forces=true)

    # Run the simulation
    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();