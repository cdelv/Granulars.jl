include("../../src/Granulars.jl")

function main()
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    g = [0.0,-9.8,0.0]
    walls = Create_Box(Lx,Ly,Lz)

    p = Particle(r=[3,1,10], v=[0,0,0], w=[0,0,-5])
    particles = [p]

    # Estimate a good time step
    dt = 0.7*PWaveTimeStep(particles)
    vis_steps = 700
    frames = 200
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, walls=walls, en=0.8, g=g)

    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();