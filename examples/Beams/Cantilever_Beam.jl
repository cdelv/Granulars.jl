include("../../src/Granulars.jl")

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]
    
    # Number of particles
    n = 20
    rad = 0.5
    m = 0.1

    # Spacing between particles
    dx = 2*rad - 0.6*rad
    
    # Create config
    walls = walls = Create_Box(25,25,25)
    conf = Config(t, dt, g=g, walls=walls, E=800000, G=1200000)
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 10, 10], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m)
        push!(particles, p)
    end
    
    # Run the simulation
    Propagate(particles, conf, vis_steps=500, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1], static=true)
end

@time main(200*500*0.0001);