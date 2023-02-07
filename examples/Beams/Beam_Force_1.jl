include("../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-0.0,0.0]
    
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # Number of particles in each direction
    # N = nx*ny*nz
    nx = 100

    # Spacing between particles
    L = 10
    dx = L/(2*nx) 
    
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
    
    # Create config
    conf = Config(t, dt, walls=walls, g=g)

    aos = Particle[]

    for i in 1:nx
        push!(aos, Particle(r=[5+i*dx, Ly/2, Lz/2], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=dx/1.5))
    end

    p1 = Particle(r=[12,Ly/2,Lz/2], v=[0,0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13,Ly/2+0.5,Lz/2], v=[0.0,0,0], w=[0,0,0], rad=2.0)
    
    # Run the simulation
    Propagate(aos, conf, vis_steps=1, file="Paraview/data", save=true, beam_forces=true)
end

@time main(200*1*0.0001);