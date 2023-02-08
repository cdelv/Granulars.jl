
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

    q0::Quaternion{Float64}
    q0ij::Quaternion{Float64}
end

"""
Convenience constructor for Beam 
- p_i:
- p_j:
- conf:

Formula for the corss section rad:
https://mathworld.wolfram.com/Sphere-SphereIntersection.html

quaternion reference frame difference:
https://math.stackexchange.com/questions/1884215/how-to-calculate-relative-pitch-roll-and-yaw-given-absolutes
"""
function Beam(p_i::Particle, p_j::Particle, conf::Config)::Beam
    L::Float64 = norm(p_i.r - p_j.r)
    r::Float64 = sqrt((-L + p_i.rad - p_j.rad)*(-L - p_i.rad + p_j.rad)*(-L + p_i.rad + p_j.rad)*(L + p_i.rad + p_j.rad))/(2*L)
    q0 = p_i.q * (-conj(Beam_Orientation(unitary(p_j.r - p_i.r)))) # goes from i to j
    qij0 = p_j.q * (-conj(p_i.q))
    Beam(conf.E, conf.G, π*r*r, L, q0, qij0)
end


function Beam_Orientation(n::SVector{3, Float64})::Quaternion{Float64}
    ex::SVector{3} = SVector(1.0,0.0,0.0) # x axis
    Angle::Float64 = angle(ex,n) # rotation angle
    u::SVector{3} = cross(ex,n)/(sin(Angle)+1e-9) # rotation axis
    angleaxis_to_quat(EulerAngleAxis(Angle, u))
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

function sparcify(v::Float64)::Float64
    if abs(v)<1e-11
        return 0.0
    else 
        return v
    end
end

"""
DESCRIPTION:

- particles: StructArray of particles.
- i:
- j:
- k:
- n:
- conf: Simulation configuration, it's a Conf struct, implemented in Configuration.jl. 

quaternion reference frame difference
https://math.stackexchange.com/questions/1884215/how-to-calculate-relative-pitch-roll-and-yaw-given-absolutes
"""
function Beam_Force(particles::StructVector{Particle}, 
    beams::StructVector{Beam}, i::Int64, j::Int64, k::Int64, 
    n::SVector{3,Float64}, conf::Config)

    # Beam Orientation 
    qb = -conj(beams.q0[k]) * particles.q[i]

    # Vector that goes from the i particle to the j particle in the Beam frame
    # Substracting (L,0,0) gets the displacement diference
    Δr = Lab_to_body(particles.r[j]-particles.r[i], qb) - SVector(beams.L[k], 0.0, 0.0)

    # Angle displacement
    Δq = particles.q[j] * (-conj(particles.q[i]))
    Δϕ = quat_to_angle(-conj(beams.q0ij[k]) * Δq, :XYZ)
    #println(quat_to_angle(beams.q0ij[k], :XYZ))
    #println(quat_to_angle(Δq, :XYZ))
    #println(Δϕ)
    #println("")

    # Create the 12x12 transformation matrix and calculate forces and torques.
    @inbounds Δs::SVector{12,Float64} = 0.5*SVector(Δr[1], Δr[2], Δr[3], Δϕ.a1, Δϕ.a2, Δϕ.a3, -Δr[1], -Δr[2], -Δr[3], -Δϕ.a1, -Δϕ.a2, -Δϕ.a3)
    Δs = sparcify.(Δs)
    F::SVector{12,Float64} = K_beam(beams.A[k], beams.L[k], conf.E, conf.G)*Δs #cte k matrix, maybe compute once

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Body_to_lab(SVector(F[1], F[2], F[3]), qb)/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Body_to_lab(SVector(F[4], F[5], F[6]), qb), particles.q[i])

    @inbounds particles.a[j] += Body_to_lab(SVector(F[7], F[8], F[9]), qb)/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Body_to_lab(SVector(F[10],F[11],F[12]), qb), particles.q[j])

    #=
    γ = 0.1
    @inbounds particles.a[i] += -γ*particles.v[i]
    @inbounds particles.τ[i] += -γ*particles.w[i]

    @inbounds particles.a[j] += -γ*particles.v[j]
    @inbounds particles.τ[j] += -γ*particles.w[j]
    =#
    return nothing
end