include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    
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
    conf = Config(t, dt, g=g, walls=walls, en=0.8, mu=0.4)

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p3 = Particle(r=[11.5, 13, 10], v=[10.0,0.0,0], w=[0,0,0], rad=2.0)

    # Run the simulation
    Propagate([p1,p2,p3], conf, vis_steps=400, file="Paraview/data", 
        save=true, beam_forces=true)
end

@time main(200*400*0.0001);