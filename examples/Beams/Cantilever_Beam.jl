include("../../src/Granulars.jl")

function main()
    # Simulation parameters
    g = [0.0,-10.0,0.0]
    
    # Number of particles
    n = 50
    rad = 0.5
    m = 0.01

    # Spacing between particles
    dx = 2*rad - 0.6*rad
    
    # Create config
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 40, 5], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m, E=800000, ν=-0.5)
        push!(particles, p)
    end

    # Estimate a good time step
    vis_steps = 100000
    frames = 200
    dt = 0.00005
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, beam_forces=true, beam_damping=true, ζ=0.1)
    
    # Run the simulation
    particles = Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", 
        save=true, fixed_spheres=[1], static=true)

    println("X,Y")
    for i in eachindex(particles)
        println(particles.r[i][1],",",particles.r[i][2])
    end
end

@time main();