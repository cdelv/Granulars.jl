include("../../src/Granulars.jl")

"""
###########################################
# REMEMBER TO UNCOMENT DAMPING IN BEAMS.JL
###########################################
"""

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,-10.0,0.0]
    
    # Number of particles
    n = 50
    rad = 0.5
    m = 0.01

    # Spacing between particles
    dx = 2*rad - 0.6*rad
    
    # Create config
    walls = walls = Create_Box(50,50,50)
    conf = Config(t, dt, g=g, walls=walls)
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 45, 5], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m, E=800000, ν=0.2)
        push!(particles, p)
    end
    
    # Run the simulation
    particles = Propagate(particles, conf, vis_steps=3000, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1], static=true)

    println("X,Y")
    for i in eachindex(particles)
        println(particles.r[i][1],",",particles.r[i][2])
    end
end

@time main(200*3000*0.0001);