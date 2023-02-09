include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.000001
    g = [0.0,0.0,0.0]
    
    # Create config
    conf = Config(t, dt, g=g, E=100000, G=100000)

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,-15.0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p3 = Particle(r=[16, 10, 10], v=[0.0,15.0,0], w=[0,0,0], rad=2.0)

    # Run the simulation
    Propagate([p1,p2,p3], conf, vis_steps=20000, file="Paraview/data", 
        save=true, beam_forces=true)
end

@time main(200*20000*0.000001);