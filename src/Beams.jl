
"""
Structure that stores information of the beams. The beams are assumed cilindrical. 
- E: Young Modulus.
- G: Shear Modulus. 
- J: Second area moment.
- I: Torsion constant.
- A: Cross section.
- L: Length of the beam.
- Δq0i: 
- Δq0j: 
"""
struct Beam
    E::Float64
    G::Float64
    J::Float64
    I::Float64 # Cilindrical Beam, then Iy = Iz.

    A::Float64
    L::Float64

    Δq0i::Quaternion{Float64}
    Δq0j::Quaternion{Float64}
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

    # we assume a cilindrical beam
    A::Float64 = π*r*r
    J::Float64 = A*A/(2*π) # π*r^4/2 -> Torsion constant.
    I::Float64 = J/2      # π*r^4/4 -> Second area moment. Iy = Iz.
    
    # Diference between the particle and beam frame. They are fixed and move the same amount.
    Δq0i::Quaternion{Float64} = Beam_Orientation(unitary(p_j.r - p_i.r)) ∘ inv(p_i.q) # Defined in Utils.jl

    # Diference between the particle and beam frame. They are fixed and move the same amount.
    Δq0j::Quaternion{Float64} = Beam_Orientation(unitary(p_j.r - p_i.r)) ∘ inv(p_j.q) # Defined in Utils.jl

    Beam(conf.E, conf.G, J, I, A, L, unitary(Δq0i), unitary(Δq0j))
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

    k::Int64 = 1

    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds j::Int64 = max(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds d::Float64 = pair[3]

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
DESCRIPTION:

- particles: StructArray of particles.
- i:
- j:
- k:

quaternion reference frame difference
https://math.stackexchange.com/questions/1884215/how-to-calculate-relative-pitch-roll-and-yaw-given-absolutes
"""
function Beam_Force(particles::StructVector{Particle}, 
    beams::StructVector{Beam}, i::Int64, j::Int64, k::Int64)
    
    # Beam Orientation in the i particle
    @inbounds qbi::Quaternion{Float64} = unitary(beams.Δq0i[k] ∘ particles.q[i])

    # Beam Orientation in the j particle
    @inbounds qbj::Quaternion{Float64} = unitary(beams.Δq0j[k] ∘ particles.q[j])

    # Vector that goes from the i particle to the j particle in the Beam frame
    # Substracting (L,0,0) gets the displacement diference
    @inbounds Δri::SVector{3, Float64} = Lab_to_body(particles.r[j]-particles.r[i], qbi) - SVector(beams.L[k], 0.0, 0.0)
    @inbounds Δrj::SVector{3, Float64} = Lab_to_body(particles.r[j]-particles.r[i], qbj) - SVector(beams.L[k], 0.0, 0.0)
    Δr::SVector{3, Float64} = 0.5*(Δri+Δrj)
    #Δr::SVector{3, Float64} = 0*Δri

    # Angle displacement
    Δq::Quaternion{Float64} = unitary(qbj ∘ inv(qbi))
    Δϕ::EulerAngles{Float64} = quat_to_angle(Δq, :XYZ)  
    #Δϕ::EulerAngles{Float64} = EulerAngles(0.0,0.0,0.0, :XYZ)     

    # Create the 12x12 transformation matrix and calculate forces and torques.
    #@inbounds Δs::SVector{12,Float64} = 0.5*SVector(Δr[1], Δr[2], Δr[3], Δϕ.a1, Δϕ.a2, Δϕ.a3, -Δr[1], -Δr[2], -Δr[3], -Δϕ.a1, -Δϕ.a2, -Δϕ.a3)
    #F::SVector{12,Float64} = K_beam(beams.A[k], beams.L[k], conf.E, conf.G)*Δs #cte k matrix, maybe compute once

    # Analitic sol of displacement times the stifness matrix.
    @inbounds F::SVector{12,Float64} = (beams.E[k]/(beams.L[k]*beams.L[k]*beams.L[k]))*SVector(
        beams.A[k]*Δr[1]*beams.L[k]*beams.L[k],
        12*Δr[2]*beams.I[k],
        12*Δr[3]*beams.I[k],
        Δϕ.a1*beams.G[k]*beams.J[k]*beams.L[k]*beams.L[k]/beams.E[k],
        Δϕ.a2*beams.I[k]*beams.L[k]*beams.L[k] - 6*Δr[3]*beams.I[k]*beams.L[k],
        Δϕ.a3*beams.I[k]*beams.L[k]*beams.L[k] + 6*Δr[2]*beams.I[k]*beams.L[k],
        -beams.A[k]*Δr[1]*beams.L[k]*beams.L[k],
        -12*Δr[2]*beams.I[k],
        -12*Δr[3]*beams.I[k],
        -Δϕ.a1*beams.G[k]*beams.J[k]*beams.L[k]*beams.L[k]/beams.E[k],
        -Δϕ.a2*beams.I[k]*beams.L[k]*beams.L[k] - 6*Δr[3]*beams.I[k]*beams.L[k],
        -Δϕ.a3*beams.I[k]*beams.L[k]*beams.L[k] + 6*Δr[2]*beams.I[k]*beams.L[k],
        )

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Body_to_lab(SVector(F[1], F[2], F[3]), qbi)/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Body_to_lab(SVector(F[4], F[5], F[6]), qbi), particles.q[i])

    @inbounds particles.a[j] += Body_to_lab(SVector(F[7], F[8], F[9]), qbi)/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Body_to_lab(SVector(F[10],F[11],F[12]), qbi), particles.q[j])

    #=
    γ = 0.5
    @inbounds particles.a[i] += -γ*particles.v[i]
    @inbounds particles.τ[i] += -γ*particles.w[i]

    @inbounds particles.a[j] += -γ*particles.v[j]
    @inbounds particles.τ[j] += -γ*particles.w[j]
    =#
    
    return nothing
end


"""
NOT USED

Computes the stifness matrix of the beam element
- A: Cross section.
- L: Length of the beam.
- E: Young Modulus.
- G: Shear Modulus. 

THIS IS NOT USED ANY MORE!
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