include("../../../src/Granulars.jl")
using Distributions

function Actions_After_Time_Step(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    cundall_particles::ExtendableSparseMatrix{Float64, Int64},
    cundall_walls::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam},
    fixed_spheres::Vector{Int64},
    static::Bool,
    t::Float64)

    nothing
end

function main(t)
    # Simulation parameters
    dt = 0.00001
    g = [0.0,-60.0,0.0]

    # Box dimensions
    Lx = 6.2
    Ly = 25.0
    Lz = 6.2
    walls = Create_Box(Lx, Ly, Lz, E=1e8, G=5e7)
    
    # Create config
    conf = Config(t, dt, g=g, walls=walls, en=0.6, mu=0.4)

    # Start particles in organized way
    particles = Particle[]
    nx = 3
    ny = 11
    nz = 3

    #=
    Lx = 8.2
    Ly = 8.2
    nx = 4
    ny = 5
    nz = 4
    =#

    dx = 2.01
    dy = 2.01
    dz = 2.01
    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz-1 
                rr = [1.1 + i*dx, 2.0 + j*dy, 1.1 + k*dz]
                push!(particles, Particle(r=rr, v=3*rand(3), rad=rand(Uniform(0.5, 1.0)), E=5e3, G=4e3))
                #push!(particles, Particle(r=rr, v=3*rand(3), rad=1.0, E=5e3, G=4e3))
            end
        end
    end

    #CSV file header
    println("y,F")

    # Run the simulation
    Propagate(particles, conf, vis_steps=1800, file="Paraview/data", 
        save=true, beam_forces=false)
end

main(200*1800*0.00001);