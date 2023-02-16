include("../../src/Granulars.jl")

function main(t)
    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25

    # time step
    dt = 0.001

    walls = Create_Box(Lx,Ly,Lz)
    walls = [Wall([0.0,1.0,0.0], [0.0,0.0,0.0])]

    conf = Config(t, dt, walls=walls)

    p = Particle(r=[10,10,10])

    Propagate([p], conf, vis_steps=200, file="Paraview/data", save=true)
end

@time main(200*0.001*200);