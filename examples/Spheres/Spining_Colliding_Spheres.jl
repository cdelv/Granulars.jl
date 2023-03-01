include("../../src/Granulars.jl")

function main()
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    g = [0.0,0.0,0.0]
    walls = Create_Box(Lx,Ly,Lz)

    p1 = Particle(r=[7,10,10], v=[5,0,0], w=[0,0,0])
    p2 = Particle(r=[18,10,10], v=[-5,0,0], w=[0,0,10])
    particles = [p1,p2]

    # Estimate a good time step
    dt = 0.7*PWaveTimeStep(particles)
    vis_steps = 110
    frames = 200
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, walls=walls, g=g, en=0.7)


    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();