include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]
    
    # Create config
    conf = Config(t, dt, g=g, gamma=10)

    p1 = Particle(r=[10, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[11, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)

    # Run the simulation
    Propagate(aos, conf, vis_steps=500, file="Paraview/data", save=true, beam_forces=true)
end

@time main(200*500*0.0001);