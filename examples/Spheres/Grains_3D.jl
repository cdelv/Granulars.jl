include("../src/Granulars.jl")

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
    
    # X coordinate walls
    W1 = Wall([1,0,0],[0,0,0])
    W2 = Wall([-1,0,0],[Lx,0,0])
    
    # Y coordinate walls
    W3 = Wall([0,1,0],[0,0,0])
    W4 = Wall([0,-1,0],[0,Ly,0])
    
    # Z coordinate walls
    W5 = Wall([0,0,1],[0,0,0])
    W6 = Wall([0,0,-1],[0,0,Lz])
    
    walls = [W1,W2,W3,W4,W5,W6]

    aos = Particle[]
    
    # Start particles in organized way
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz 
                push!(aos, Particle([i*dx, j*dy, k*dz], rand(3)))
            end
        end
    end
    
    # Create config
    conf = Config(t, dt, walls=walls, g=g, gamma=5)
    
    # Run the simulation
    Propagate(aos, conf, vis_steps=60, file="Paraview/data", save=true)
end

@time main(200*60*0.001);