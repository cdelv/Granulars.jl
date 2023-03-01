include("../../src/Granulars.jl")

function main()
    # This simulates a Spinning T, see video:
    # https://www.youtube.com/watch?v=1n-HMSCDYtM

    # Simulation parameters
    g = [0.0,0.0,0.0]
    
    # Spinning speed
    v = 66.6

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,-v,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13.8, 10, 10], v=[0.0,0.0,0], w=[0,0.0,v/4], rad=2.0)
    p3 = Particle(r=[17.6, 10, 10], v=[0.0,v,0], w=[0,0,0], rad=2.0)
    p4 = Particle(r=[13.8, 10, 13.8], v=[0.0,0.01,0], w=[0,0.0,v/4], rad=2.0)
    particles = [p1,p2,p3,p4]

    # Estimate a good time step
    vis_steps = 600
    frames = 200
    dt = 0.00005
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, beam_forces=true)

    # Run the simulation
    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();