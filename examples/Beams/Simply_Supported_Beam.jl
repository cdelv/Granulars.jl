include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]
    
    # Number of particles
    n = 21
    rad = 0.5
    m = 0.1

    # Spacing between particles
    dx = 2*rad - 0.3*rad
    
    # Create config
    conf = Config(t, dt, g=g, E=80000, G=120000)
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 10, 10], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m)
        push!(particles, p)
    end
    
    # Run the simulation
    Propagate(particles, conf, vis_steps=500, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1,n], static=true)
end

@time main(200*500*0.0001);