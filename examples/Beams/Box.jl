include("../../src/Granulars.jl")

function main()
    # Simulation parameters
    g = [0.0,-9.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx,Ly,Lz, E=1e9, ν=-0.5)
    
    q = angle_to_quat(EulerAngles(0.5,0.5,0.5, :XYZ))

    # Start particles in organized way
    particles = Particle[]
    nx = 4
    ny = 4
    nz = 4
    dx = 1.75
    dy = 1.75
    dz = 1.75
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz 
                rr = Lab_to_body(SVector(5.0 + i*dx, 5.5 + j*dy, 9.0 + k*dz), q)
                push!(particles, Particle(r=rr, E=1e4, ν=0.2, v=[1.0,0.0,0.0]))
            end
        end
    end

    # Estimate a good time step
    dt = 0.4*PWaveTimeStep(particles)
    vis_steps = 30
    frames = 200
    println("dt = ", dt)

    # Create config
    t = frames*vis_steps*dt
    conf = Config(t, dt, g=g, walls=walls, mu=0.4, en=0.85, v=0.01, beam_forces=true, 
        thorsten_damping=true, beam_damping=true, ζ=0.05, fracture=true, c=1300.0)

    # Run the simulation
    Propagate(particles, conf, vis_steps=vis_steps, file="Paraview/data", save=true)
end

@time main();