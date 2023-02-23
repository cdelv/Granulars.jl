include("../../src/Granulars.jl")

function main(t)
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # time step
    dt = 0.001

    walls = Create_Box(Lx,Ly,Lz)

    conf = Config(t, dt, walls=walls, en=0.9, v=1.0, mu=0.6)

    p1 = Particle(r=[10,10,10], w=[0,0,0])

    Propagate([p1], conf, vis_steps=200, file="Paraview/data", save=true)
end

@time main(200*0.001*200);