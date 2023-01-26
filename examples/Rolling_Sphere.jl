include("../src/Granulars.jl")

function main(t)
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # time step
    dt = 0.01

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

    conf = Config(100, dt, walls=walls)

    p = Particle(r=[3,1,10], v=[5,0,0], w=[0,0,0])

    Propagate([p], conf, vis_steps=65, file="Paraview/data", save=true)
end

@time main(100);