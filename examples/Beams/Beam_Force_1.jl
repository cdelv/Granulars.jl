include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]
    
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # Number of particles in each direction
    # N = nx*ny*nz
    nx = 20

    # Spacing between particles
    L = 15
    dx = 1.1*L/(nx) 
    
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
    conf = Config(t, dt, walls=walls, g=g, gamma=10)

    aos = Particle[]

    for i in 1:nx
        p = Particle(r=[5+i*dx, Ly/2, Lz/2], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=0.8)
        p = Set_q(p,Quaternion(1.0I))
        push!(aos, p)
    end

    push!(aos, Particle(r=[5+10*dx, Ly/2+8, Lz/2], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=0.8))

    p1 = Particle(r=[12   ,Ly/2    ,Lz/2], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[15.5 ,Ly/2+1    ,Lz/2], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p3 = Particle(r=[19   ,Ly/2+2  ,Lz/2], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    
    # Run the simulation
    Propagate(aos, conf, vis_steps=500, file="Paraview/data", save=true, beam_forces=true)
end

@time main(200*500*0.0001);