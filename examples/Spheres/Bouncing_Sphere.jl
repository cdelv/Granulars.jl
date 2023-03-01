include("../../src/Granulars.jl")

function main()
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx,Ly,Lz)
    g = [0, -9.8, 0]

    p1 = Particle(r=[10,10,10], w=[0,0,0])
    particles=[p1]

    # Estimate a good time step
    dt = 0.7*PWaveTimeStep(particles)
    vis_steps = 200
    frames = 200
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, walls=walls, g=g, en=0.8, v=0.8, mu=0.6, thorsten_damping=true)

    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();