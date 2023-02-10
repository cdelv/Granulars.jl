include("../../src/Granulars.jl")

function main(t)
    # This simulates a Spinning T, see video:
    # https://www.youtube.com/watch?v=1n-HMSCDYtM

    # Simulation parameters
    dt = 0.0001
    g = [0.0,0.0,0.0]
    
    # Create config
    conf = Config(t, dt, g=g, E=500000, G=500000)

    # Spinning speed
    v = 66.6

    # Create particles
    p1 = Particle(r=[10, 10, 10], v=[0.0,-v,0], w=[0,0,0], rad=2.0)
    p2 = Particle(r=[13.8, 10, 10], v=[0.0,0.0,0], w=[0,0.0,v/4], rad=2.0)
    p3 = Particle(r=[17.6, 10, 10], v=[0.0,v,0], w=[0,0,0], rad=2.0)
    p4 = Particle(r=[13.8, 10, 13.8], v=[0.0,0.0,0], w=[0,0.0,v/4], rad=2.0)

    # Run the simulation
    Propagate([p1,p2,p3,p4], conf, vis_steps=400, file="Paraview/data", 
        save=true, beam_forces=true)
end

@time main(200*400*0.0001);