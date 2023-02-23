include("../../src/Granulars.jl")
l = 1.0
g = 1.0

function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, conf::Config,
    cundall_particles::ExtendableSparseMatrix{Float64, Int64},
    cundall_walls::ExtendableSparseMatrix{Float64, Int64},
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::StructVector{Beam})
    q = particles.q[1]
    
    particles.τ[1] = particles.m[1]*g*l*SVector(
            2.0*(q.q2*q.q3 + q.q0*q.q1),
            -2.0*(q.q1*q.q3 - q.q0*q.q2), 
            0.0)
    
    return nothing
end

function main(t)
    p = Particle(w=[0, 1.0, 9.0])
    p = Set_I(p,SVector{3,Float64}([0.2, 0.2, 1.0]))
    p = Set_q(p,angle_to_quat(π/7, 0, -π/2,:XYZ))
    conf = Config(t,0.0001)
    Propagate([p], conf; vis_steps=20, file="Paraview/data", save=true)
end

main(0.5)