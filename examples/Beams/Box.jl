include("../../src/Granulars.jl")

dt = 0.00005

function main(t)
    # Simulation parameters
    g = [0.0,-9.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx,Ly,Lz, E=1e9, ν=-0.5)
    
    # Create config
    conf = Config(t, dt, g=g, walls=walls, en=1.0, v=0.3, thorsten_damping=false, beam_damping=true)

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
                rr = Lab_to_body(SVector(9.0 + i*dx, 4.5 + j*dy, 9.0 + k*dz), q)
                push!(particles, Particle(r=rr, E=9e3, ν=-0.2))
            end
        end
    end

    global dt = 0.7*PWaveTimeStep(particles)
    println("dt = ", dt)

    # Run the simulation
    Propagate(particles, conf, vis_steps=1600, file="Paraview/data", save=true, beam_forces=true)
end

@time main(200*1600*dt);