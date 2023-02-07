
"""
Structure that stores information of the beams
- E: Young Modulus.
- G: Shear Modulus. 
- A: Cross section.
- L: Length of the beam.
- r_i: 
- r_j: 
- q_i: 
- q_j: 

TO DO: MAKE IT WORK
"""
struct Beam
    E::Float64
    G::Float64

    A::Float64
    L::Float64

    r_i::SVector{3, Float64}
    r_j::SVector{3, Float64}

    q_i::Quaternion{Float64}
    q_j::Quaternion{Float64}
end

"""
Convenience constructor for Beam 
- p_i:
- p_j:
- conf:

Formula for the corss section rad:
https://mathworld.wolfram.com/Sphere-SphereIntersection.html
"""
function Beam(p_i::Particle, p_j::Particle, conf::Config)::Beam
    L::Float64 = norm(p_i.r - p_j.r)
    r::Float64 = sqrt((-L + p_i.rad - p_j.rad)*(-L - p_i.rad + p_j.rad)*(-L + p_i.rad + p_j.rad)*(L + p_i.rad + p_j.rad))/(2*L)
    Beam(conf.E, conf.G, L, π*r*r, p_i.r, p_j.r, p_i.q, p_j.q)
end

"""
Creates beams between all intersecting particles
- neighborlist:
- conf:
- beam_bonds:
- beams:
"""
function Create_beams(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    conf::Config,
    beam_bonds::ExtendableSparseMatrix{Int64, Int64},
    beams::Vector{Beam})

    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds j::Int64 = max(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds d::Float64 = pair[3]
        k::Int64 = 1

        # Interpenetration distance.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s > 0.0
            beam_bonds[i,j] = k
            push!(beams, Beam(particles[i], particles[j], conf))
            k+=1
        end
    end

    return nothing
end

"""
Computes the stifness matrix of the beam element
- A: Cross section.
- L: Length of the beam.
- E: Young Modulus.
- G: Shear Modulus. 
"""
function K_beam(A::Float64, L::Float64, E::Float64, G::Float64)::SMatrix{12, 12, Float64, 144}
    # we assume a cilindrical beam
    J::Float64 = A*A/(2*π) # π*r^4/2 -> Torsion constant.
    Iy::Float64 = J/2      # π*r^4/4 -> Second area moment.
    Iz::Float64 = Iy       # π*r^4/4 -> Second area moment.
    return SMatrix{12,12}(
        A*L*L, 0.0, 0.0, 0.0, 0.0, 0.0, -A*L*L, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 12*Iz, 0.0, 0.0, 0.0, 6*L*Iz, 0.0, -12*Iz, 0.0, 0.0, 0.0, 6*L*Iz,
        0.0, 0.0, 12*Iy, 0.0, -6*L*Iy, 0.0, 0.0, 0.0, -12*Iy, 0.0, -6*L*Iy, 0.0,
        0.0, 0.0, 0.0, G*J*L*L/E, 0.0, 0.0, 0.0, 0.0, 0.0, -G*J*L*L/E, 0.0, 0.0,
        0.0, 0.0, -6*L*Iy, 0.0, 4*L*L*Iy, 0.0, 0.0, 0.0, 6*L*Iy, 0.0, 2*L*L*Iy, 0.0,
        0.0, 6*L*Iz, 0.0, 0.0, 0.0, 4*L*L*Iz, 0.0, -6*L*Iz, 0.0, 0.0, 0.0, 2*L*L*Iz,
        -A*L*L, 0.0, 0.0, 0.0, 0.0, 0.0, A*L*L, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, -12*Iz, 0.0, 0.0, 0.0, -6*L*Iz, 0.0, 12*Iz, 0.0, 0.0, 0.0, -6*L*Iz,
        0.0, 0.0, -12*Iy, 0.0, 6*L*Iy, 0.0, 0.0, 0.0, 12*Iy, 0.0, 6*L*Iy, 0.0,
        0.0, 0.0, 0.0, -G*J*L*L/E, 0.0, 0.0, 0.0, 0.0, 0.0, G*J*L*L/E, 0.0, 0.0,
        0.0, 0.0, -6*L*Iy, 0.0, 2*L*L*Iy, 0.0, 0.0, 0.0, 6*L*Iy, 0.0, 4*L*L*Iy, 0.0,
        0.0, 6*L*Iz, 0.0, 0.0, 0.0, 2*L*L*Iz, 0.0, -6*L*Iz, 0.0, 0.0, 0.0, 4*L*L*Iz
    )*(E/L^3)
end


"""
DESCRIPTION:

- particles: StructArray of particles.
- i:
- j:
- k:
- n:
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 
"""
function Beam_Force(particles::StructVector{Particle}, 
    beams::StructVector{Beam}, i::Int64, j::Int64, k::Int64, 
    n::SVector{3,Float64}, conf::Config)

    # Calculate displacements in the beam frame.
    dxi::SVector{3,Float64} = Lab_to_Beam(n, beams.r_i[k] - particles.r[i])
    dΩi::SVector{3,Float64} = Lab_to_Beam(n, quat_to_angle(beams.q_i[k]) - quat_to_angle(particles.q[i]))
    dxj::SVector{3,Float64} = Lab_to_Beam(n, beams.r_j[k] - particles.r[j])
    dΩj::SVector{3,Float64} = Lab_to_Beam(n, quat_to_angle(beams.q_j[k]) - quat_to_angle(particles.q[j]))

    # Create the 12x12 transformation matrix and calculate forces and torques.
    @inbounds Δs::SVector{12,Float64} = SVector(dxi[1], dxi[2], dxi[3], dΩi[1], dΩi[2], dΩi[3], dxj[1], dxj[2], dxj[3], dΩj[1], dΩj[2], dΩj[3])
    F::SVector{12,Float64} = K_beam(beams.A[k], beams.L[k], conf.E, conf.G)*Δs #cte k matrix, maybe compute once

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Beam_to_Lab(n, SVector(F[1], F[2], F[3]))/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Beam_to_Lab(n, SVector(F[4], F[5], F[6])), particles.q[i])

    @inbounds particles.a[j] += Beam_to_Lab(n, SVector(F[7], F[8], F[9]))/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Beam_to_Lab(n, SVector(F[10],F[11],F[12])), particles.q[j])

    γ = 0.1
    @inbounds particles.a[i] += -γ*particles.v[i]
    @inbounds particles.τ[i] += -γ*particles.w[i]

    @inbounds particles.a[j] += -γ*particles.v[j]
    @inbounds particles.τ[j] += -γ*particles.w[j]

    # Update the beam information (NOT NEEDED DUE TO THE WAY THE K MATRIX WORKS)
    #beams.r_i[k] = particles.r[i]
    #beams.r_j[k] = particles.r[j]
    #beams.q_i[k] = particles.q[i]
    #beams.q_j[k] = particles.q[j]

    return nothing
end