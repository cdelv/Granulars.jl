include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.00005
    g = [0.0,-9.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx,Ly,Lz)
    
    # Create config
    conf = Config(t, dt, g=g, walls=walls, E=10000, G=10000, gamma=15)

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
                push!(particles, Particle(r=[rr[1],rr[2],rr[3]]))
            end
        end
    end

    # Run the simulation
    Propagate(particles, conf, vis_steps=1600, file="Paraview/data", 
        save=true, beam_forces=true)
end

@time main(200*1600*0.00005);