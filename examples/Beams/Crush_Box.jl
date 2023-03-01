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
    conf = Config(t, dt, g=g, walls=walls, en=0.85, v=0.01, thorsten_damping=false, beam_damping=true, ζ=0.06)

    # Start particles in organized way
    particles = Particle[]
    nx = 4*2
    ny = 3*2
    nz = 4*2
    dx = 1.75
    dy = 1.75
    dz = 1.75
    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz -1
                rr = SVector(9.0 + i*dx/2, 1.0/2 + j*dy/2, 9.0 + k*dz/2)
                push!(particles, Particle(r=rr, E=9e3, ν=-0.2, rad=0.5))
            end
        end
    end

    push!(particles, Particle(r=[9.0 + 1.5*dx, 1.0 + 5*dy, 9.0 + 1*dz], v=[0, -0.6, 0], E=1e5, ν=0.2))

    global dt = 0.7*PWaveTimeStep(particles)
    println("dt = ", dt)

    # Run the simulation
    Propagate(particles, conf, vis_steps=1600, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[length(particles)])
end

@time main(200*1600*dt);