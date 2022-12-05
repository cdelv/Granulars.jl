include("../src/Granulars.jl")
function main()
    t = 30
    dt = 0.01
    g = [0,-1,0]
    
    Lx = 30
    Ly = 30
    Lz = 30
    
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
    
    N = 20
    aos = [Particle(14*rand(3),rand(3)) for i = 1:N]
    
    conf = Config(t,dt,g,walls=walls)
    
    Propagate(aos, conf, vis_steps=50, file="Paraview/data", save=true)
end

@time main()
@time main();