include("../../../../src/Granulars.jl")

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
    
    v::Float64 = 1.0e-3
    y::Float64 = 8.9 - v*t
    conf.walls[4] = Set_Q(conf.walls[4], SVector(0.0, y ,0.0))

    println(y,",", conf.walls[4].F[2])

    nothing
end

function main(t)
    # Simulation parameters
    dt = 0.0001
    g = [0.0,0.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 25
    Lz = 25
    walls = Create_Box(Lx, Ly, Lz, E=1e15, G=1e15)
    
    # Create config
    conf = Config(t, dt, g=g, walls=walls, en=0.8)

    # Start particles in organized way
    particles = Particle[]
    nx = 4
    ny = 5
    nz = 4
    dx = 1.75
    dy = 1.75
    dz = 1.75
    for i in 0:nx-1
        for j in 0:ny-1
            for k in 0:nz-1 
                rr = [9.0 + i*dx, 1.0 + j*dy, 9.0 + k*dz]
                push!(particles, Particle(r=rr, E=1e7, G=1e7))
            end
        end
    end

    #CSV file header
    println("y,F")

    # Run the simulation
    Propagate(particles, conf, vis_steps=5000, file="Paraview/data", 
        save=true, beam_forces=true)
end

main(200*5000*0.0001);