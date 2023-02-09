include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.00001
    g = [0.0,-10.0,0.0]
    
    # Number of particles
    n = 20
    rad = 0.5
    m = 0.1

    # Spacing between particles
    dx = 2*rad - 0.6*rad
    
    # Create config
    conf = Config(t, dt, g=g, E=800000, G=1200000)
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 10, 10], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m)
        push!(particles, p)
    end
    
    # Run the simulation
    Propagate(particles, conf, vis_steps=5000, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1], static=true)
end

@time main(200*5000*0.00001);