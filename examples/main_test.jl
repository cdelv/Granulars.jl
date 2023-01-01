include("../src/Granulars.jl")
function main(t,N)
    dt = 0.01
    g = [0.0,-9.0,0.0]
    
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
    
    aos = [Particle((Lx-2.0)*rand(3)+[1,1,1],4.0*rand(3)) for i = 1:N]
    
    conf = Config(t,dt,g=g,walls=walls)
    
    Propagate(aos, conf, vis_steps=65, file="Paraview/data", save=true)
end

@time main(1,20);