include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,0.0,0.0]
    
    # Create config
    conf = Config(t, dt, g=g, E=1000, G=1000)

    p1 = Particle(r=[10, 10, 10], v=[0.0,0.0,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13, 10, 10], v=[10.0,0.0,0], w=[0,0,0], rad=2.0)

    # Run the simulation
    Propagate([p1,p2], conf, vis_steps=50, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1], static=true)
end

@time main(200*50*0.0001);