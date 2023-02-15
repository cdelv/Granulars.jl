include("../../src/Granulars.jl")

function main(t)
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # time step
    dt = 0.001

    walls = Create_Box(Lx,Ly,Lz)

    conf = Config(t, dt, walls=walls, en=0.8)

    p = Particle(r=[3,1,10], v=[0,0,0], w=[0,0,-5])

    Propagate([p], conf, vis_steps=600, file="Paraview/data", save=true)
end

@time main(200*0.001*600);