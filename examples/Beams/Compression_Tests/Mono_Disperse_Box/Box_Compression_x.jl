include("../../../../src/Granulars.jl")
using CSV, DataFrames

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
    y::Float64 = 8.35 - v*t
    conf.walls[4] = Set_Q(conf.walls[4], SVector(0.0, y ,0.0))

    println(y,",", conf.walls[4].F[2])

    nothing
end

function main(t)
    # Load box
    data = CSV.read("Beams/Compression_Tests/Mono_Disperse_Box/data_200.csv", DataFrame)
    x = data[!, "x"]
    y = data[!, "y"]
    z = data[!, "z"]

    # Simulation parameters
    dt = 0.00001
    g = [0.0,0.0,0.0]

    # Box dimensions
    Lx = 25
    Ly = 8.35
    Lz = 25
    walls = Create_Box(Lx, Ly, Lz, E=1e15, G=1e15)
    
    # Create config
    conf = Config(t, dt, g=g, walls=walls, en=0.6, mu=0.4)

    # Load Particles
    q = angle_to_quat(EulerAngles(0.0,0.0, π/2, :XYZ))
    particles = Particle[]

    for i in eachindex(x)
        rr = SVector(x[i], y[i], z[i])
        rr = Lab_to_body(rr, q) + SVector(4.0, 8.25, 4.0)
        push!(particles, Particle(r=rr, v=zeros(3), w=zeros(3), E=1e7, G=1e7))
    end

    #CSV file header
    println("y,F")

    # Run the simulation
    Propagate(particles, conf, vis_steps=10000, file="Paraview/data", 
        save=true, beam_forces=true)
end

main(200*10000*0.0001);