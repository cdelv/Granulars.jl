include("../../src/Granulars.jl")
using Distributions

function main(t)
    # Simulation parameters
    dt = 0.001
    g = [0.0,-9.0,0.0]
    
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # Number of particles in each direction
    # N = nx*ny*nz
    nx = 10
    ny = 10
    nz = 10

    # Spacing between particles
    dx = Lx/(nx+1) 
    dy = Ly/(ny+1)
    dz = Lz/(nz+1)
    
    # Simulation walls
    walls = Create_Box(Lx,Ly,Lz)

    particles = Particle[]
    
    # Start particles in organized way
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz 
                v = 2.0*unitary(SVector{3}(rand(Uniform(-1, 1), 3)))
                push!(particles, Particle(r=[i*dx, j*dy, k*dz], v=[v[1],v[2],v[3]]))
            end
        end
    end
    
    # Create config
    conf = Config(t, dt, walls=walls, g=g, gamma=5)
    
    # Run the simulation
    Propagate(particles, conf, vis_steps=60, file="Paraview/data", save=true)
end

@time main(200*60*0.001);