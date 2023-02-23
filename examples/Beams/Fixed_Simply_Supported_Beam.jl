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
    n = 51
    rad = 0.5
    m = 0.01

    # Spacing between particles
    dx = 2*rad - 0.6*rad
    
    # Create config
    conf = Config(t, dt, g=g)
    particles = Particle[]

    for i in 1:n
        p = Particle(r=[5+i*dx, 5, 45], v=[0.0,0.0,0.0], w=[0.0,0.0,0.0], rad=rad, m=m, E=800000, ν=-0.5)
        push!(particles, p)
    end
    
    # Run the simulation
    particles = Propagate(particles, conf, vis_steps=7000, file="Paraview/data", 
        save=true, beam_forces=true, fixed_spheres=[1,n], static=true)

    println("X,Y")
    for i in eachindex(particles)
        println(particles.r[i][1],",",particles.r[i][2])
    end
end

@time main(200*7000*0.0001);