include("../../src/Granulars.jl")

function main(t)
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # time step
    dt = 0.0001
    g = [0.0,0.0,0.0]

    walls = Create_Box(Lx,Ly,Lz)

    conf = Config(t, dt, walls=walls, g=g, en=0.7)

    p1 = Particle(r=[7,10,10], v=[5,0,0], w=[0,0,0])
    p2 = Particle(r=[18,10,10], v=[-5,0,0], w=[0,0,10])


    Propagate([p1,p2], conf, vis_steps=600, file="Paraview/data", save=true)
end

@time main(200*0.0001*600);