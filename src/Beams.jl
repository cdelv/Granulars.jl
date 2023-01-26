
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
        k::Int64 = 0

        # Interpenetration distance.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s > 0.0
            k+=1
            beam_bonds[i,j] = k
            push!(beams, Beam(particles[i], particles[j], conf))
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

# Compute beam force