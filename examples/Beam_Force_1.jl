include("../src/Granulars.jl")

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
    nx = 10

    # Spacing between particles
    L = 10
    dx = L/nx 
    
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
    conf = Config(t, dt, walls=walls, g=g, gamma=5)

    aos = Particle[]

    for i in 1:nx
        push!(aos, Particle(r=[5+i*dx, Ly/2, Lz/2], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=0.8))
    end

    p1 = Particle(r=[12,Ly/2,Lz/2], v=[0.0,0.0,0], w=[1,0,0], rad=2.0)
    p2 = Particle(r=[16-1,Ly/2,Lz/2], v=[0.0,0.0,0], w=[1,0,0], rad=2.0)
    
    p3 = Particle(r=[20-0.6,Ly/2,Lz/2], v=[0.0,-0.0,0], w=[0,0,0], rad=2.0)

    #p1 = Set_q(p1,Quaternion(1.0I))
    #p2 = Set_q(p2,Quaternion(1.0I))
    #p3 = Set_q(p3,Quaternion(1.0I))
    
    #=
    p1 = Set_I(p1,2*(2/5)*SVector(1.0,1.0,1.0))
    p2 = Set_I(p2,2*(2/5)*SVector(1.0,1.0,1.0))
    =#
    
    # Run the simulation
    Propagate([p1,p2], conf, vis_steps=300, file="Paraview/data", save=true, beam_forces=true)
end

@time main(200*300*0.0001);