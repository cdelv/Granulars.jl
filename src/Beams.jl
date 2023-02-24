"""
Structure that stores information of the beams. The beams are assumed cilindrical. 
- E: Young Modulus.
- G: Shear Modulus. 
- J: Torsion constant.
- I: Second area moment. 
- A: Cross section.
- L: Length of the beam.
- Δq0i: Initial twisting of the ith particle with respect to the beam.
- Δq0j: Initial twisting of the jth particle with respect to the beam.
- K: Stifness matrix.
- C: Damping matrix.
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

    K::SMatrix{12, 12, Float64, 144}
    C::SMatrix{12, 12, Float64, 144}
end

"""
Convenience constructor for Beam.
- p_i: Particle on one of the ends of the beam.
- p_j: Particle on one of the ends of the beam.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 

TO DO: THINK ON A BETTER CRITERIA FOR BEAM PROPERTIES

Formula for the corss section rad:
https://mathworld.wolfram.com/Sphere-SphereIntersection.html

- quaternion reference frame difference:
https://math.stackexchange.com/questions/1884215/how-to-calculate-relative-pitch-roll-and-yaw-given-absolutes

NORMAL MODES:
- Vibrations of a Free-Free Beam by Mauro Caresta:
http://www.varg.unsw.edu.au/Assets/link%20pdfs/Beam_vibration.pdf

- Continuous Systems with Longitudinal Vibration by Tom Irvine:
https://endaq.com/pages/continuous-systems-with-longitudinal-vibration

- Torsional Vibrations in Free-Free Bar with Rectangular Cross-Section by Daniel A. Russell:
https://www.acs.psu.edu/drussell/Demos/Torsional/torsional.html
"""
function Beam(p_i::Particle, p_j::Particle, conf::Config)::Beam
    L::Float64 = norm(p_i.r - p_j.r)
    r::Float64 = sqrt(4.0*L^2*p_i.rad^2 - (L^2 - p_j.rad^2 + p_i.rad^2 )^2 )/(2.0*L)

    # we assume a cilindrical beam
    A::Float64 = π*r*r
    J::Float64 = A*A/(2.0*π) # π*r^4/2 -> Torsion constant.
    II::Float64 = J/2.0       # π*r^4/4 -> Second area moment. Iy = Iz.
    
    # Diference between the particle and beam frame. They are fixed and move the same amount.
    Δq0i::Quaternion{Float64} = Beam_Orientation(unitary(p_j.r - p_i.r)) ∘ inv(p_i.q) # Defined in Utils.jl

    # Diference between the particle and beam frame. They are fixed and move the same amount.
    Δq0j::Quaternion{Float64} = Beam_Orientation(unitary(p_j.r - p_i.r)) ∘ inv(p_j.q) # Defined in Utils.jl

    # Average Young and shear modulus
    # TO DO: THINK ON A BETTER CRITERIA
    Eij::Float64 = 0.5*(p_i.E + p_j.E)
    Gij::Float64 = 0.5*(p_i.G + p_j.G)

    # Beam density
    # TO DO: THINK ON A BETTER CRITERIA
    ρ::Float64 = 0.5*(Get_Density(p_i) + Get_Density(p_j))

    # Compute Mass ans Stifness Matrix
    K::SMatrix{12, 12, Float64, 144} = Stifness_Matrix(A, L, Eij, Gij, II, J)
    M::SMatrix{12, 12, Float64, 144} = Mass_Matrix(A, L, ρ)

    # Compute first 3 normal modes on each direction
    # This is just an estimate. The real way to do it is to solve the eigenvalue problem
    # det(K - ω²M) = 0
    # But my matrices are too ilconditioned.
    ω::SVector{9, Float64} = sort(SVector(
            sqrt(Eij*II/(ρ*A))*(4.73/L)^2,    # Transversal Free-Free beam
            sqrt(Eij*II/(ρ*A))*(7.8532/L)^2,
            sqrt(Eij*II/(ρ*A))*(10.9956/L)^2,
            1.0*sqrt(Eij/ρ)/(4.0*π*L),        # Longitudinal Free-Free beam
            2.0*sqrt(Eij/ρ)/(4.0*π*L),
            3.0*sqrt(Eij/ρ)/(4.0*π*L),
            1.0*sqrt(Gij*J/(ρ*II))/(4.0*π*L), # Torsional Free-Free beam
            2.0*sqrt(Gij*J/(ρ*II))/(4.0*π*L),
            3.0*sqrt(Gij*J/(ρ*II))/(4.0*π*L)
        ))

    # Pick first and (2-3) normal modes for computing the damping coefficents
    ω1::Float64 = ω[1]
    ω2::Float64 = ω[3]

    # Compute Rayleigh damping coefficents
    a0::Float64 = 2.0*conf.ζ*ω1*ω2/(ω1+ω2)
    a1::Float64 = 2.0*conf.ζ/(ω1+ω2)

    # Compute Rayleigh Damping Matrix.
    C::SMatrix{12, 12, Float64, 144} = a0*M + a1*K

    # Create Beam
    Beam(Eij, Gij, J, II, A, L, unitary(Δq0i), unitary(Δq0j), K, C)
end

"""
Creates beams between all intersecting particles
- neighborlist: Neighbor list for particle-to-particle interaction force calculations.
- beam_bonds: Dictionary that stores wich beam connects with each particle pair.
- beams: StructArray of beams between particles.
- conf: Simulation configuration, its a Conf struct, implemented in Configuration.jl. 
"""
function Create_beams!(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    beam_bonds::Dict{Tuple{Int64, Int64}, Int64},
    beams::Vector{Beam}, conf::Config)

    k::Int64 = 1

    for pair in neighborlist
        @inbounds i::Int64 = min(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds j::Int64 = max(pair[1], pair[2]) # For symetric acces to the beam_bonds matrix
        @inbounds d::Float64 = pair[3]

        # Interpenetration distance.
        @inbounds s::Float64 = particles.rad[i] + particles.rad[j] - d

        # Check for contact. Remember that the neighborlist hass a bigger cuttof. 
        if s > 0.0
            beam_bonds[(i,j)] = k
            push!(beams, Beam(particles[i], particles[j], conf))
            k+=1
        end
    end

    return nothing
end

"""
Calculates the beam forces

- particles: StructArray of particles.
- beam_bonds: Dictionary that stores wich beam connects with each particle pair.
- i: Index of the ith particle that conforms the beam.
- j: Index of the jth particle that conforms the beam.
- k: Index of the beam that corresponds to the pair of particles (i,j).
- conf: 

TO DO: CHEK IF UNITARY(Q) IS NECESARY

quaternion reference frame difference
https://math.stackexchange.com/questions/1884215/how-to-calculate-relative-pitch-roll-and-yaw-given-absolutes
"""
function Beam_Force!(particles::StructVector{Particle}, 
    beams::StructVector{Beam}, i::Int64, j::Int64, k::Int64, conf::Config)
    
    # Beam Orientation in the i particle
    @inbounds qbi::Quaternion{Float64} = unitary(beams.Δq0i[k] ∘ particles.q[i])

    # Beam Orientation in the j particle
    @inbounds qbj::Quaternion{Float64} = unitary(beams.Δq0j[k] ∘ particles.q[j])

    # Vector that goes from the i particle to the j particle in the Beam frame
    # Substracting (L,0,0) gets the displacement diference
    @inbounds Δri::SVector{3, Float64} = Lab_to_body(particles.r[j]-particles.r[i], qbi) - SVector(beams.L[k], 0.0, 0.0)
    @inbounds Δrj::SVector{3, Float64} = Lab_to_body(particles.r[j]-particles.r[i], qbj) - SVector(beams.L[k], 0.0, 0.0)
    Δr::SVector{3, Float64} = 0.5*(Δri+Δrj)

    # Angle displacement
    Δq::Quaternion{Float64} = unitary(qbj ∘ inv(qbi))
    Δϕ::EulerAngles{Float64} = quat_to_angle(Δq, :XYZ)      

    # Create the 12x12 transformation matrix and calculate forces and torques.
    # Each extreme hass the same displacement but negative relative to the other one.
    # The reference frame is in the center of the beam. 
    @inbounds Δs::SVector{12,Float64} = 0.5*SVector(Δr[1], Δr[2], Δr[3], Δϕ.a1, Δϕ.a2, Δϕ.a3, -Δr[1], -Δr[2], -Δr[3], -Δϕ.a1, -Δϕ.a2, -Δϕ.a3)
    @inbounds F::SVector{12,Float64} = beams.K[k]*Δs
    
    # Compute relative velocities and damping force.
    if conf.beam_damping
        @inbounds v::SVector{3, Float64} = Lab_to_body(particles.v[i], qbi) - Lab_to_body(particles.v[j], qbi)
        @inbounds w::SVector{3, Float64} = Lab_to_body(Body_to_lab(particles.w[i],particles.q[i]), qbi) - Lab_to_body(Body_to_lab(particles.w[j],particles.q[j]), qbi)
        @inbounds V::SVector{12,Float64} = SVector(v[1], v[2], v[3], w[1], w[2], w[3], -v[1], -v[2], -v[3], -w[1], -w[2], -w[3])
        @inbounds F -= beams.C[k]*V
    end

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Body_to_lab(SVector(F[1], F[2], F[3]), qbi)/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Body_to_lab(SVector(F[4], F[5], F[6]), qbi), particles.q[i])
    @inbounds particles.a[j] += Body_to_lab(SVector(F[7], F[8], F[9]), qbi)/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Body_to_lab(SVector(F[10],F[11],F[12]), qbi), particles.q[j])
    
    return nothing
end


"""
Computes the stifness matrix of a beam element.
- A: Cross section.
- L: Length of the beam.
- E: Young Modulus.
- G: Shear Modulus.
- J: Torsion constant.  

Taken from:
Dynamic Analysis of Structures by John T. Katsikadelis: eq (11.9.9)
"""
function Stifness_Matrix(A::Float64, L::Float64, E::Float64, G::Float64, I::Float64, J::Float64)::SMatrix{12, 12, Float64, 144}
    # we assume a cilindrical beam
    Iy::Float64 = I      # π*r^4/4 -> Second area moment.
    Iz::Float64 = I       # π*r^4/4 -> Second area moment.
    return SMatrix{12, 12, Float64, 144}(
        L*L*A ,      0.0,       0.0,        0.0,        0.0,        0.0, -L*L*A,       0.0,      0.0,        0.0,        0.0,        0.0,
        0.0   ,  12.0*Iz,       0.0,        0.0,        0.0,   6.0*L*Iz,    0.0,  -12.0*Iz,      0.0,        0.0,        0.0,   6.0*L*Iz,
        0.0   ,      0.0,   12.0*Iy,        0.0,  -6.0*L*Iy,        0.0,    0.0,       0.0, -12.0*Iy,        0.0,  -6.0*L*Iy,        0.0,
        0.0   ,      0.0,       0.0,  G*L*L*J/E,        0.0,        0.0,    0.0,       0.0,      0.0, -G*L*L*J/E,        0.0,        0.0,
        0.0   ,      0.0, -6.0*L*Iy,        0.0, 4.0*L*L*Iy,        0.0,    0.0,       0.0, 6.0*L*Iy,        0.0, 2.0*L*L*Iy,        0.0,
        0.0   , 6.0*L*Iz,       0.0,        0.0,        0.0, 4.0*L*L*Iz,    0.0, -6.0*L*Iz,      0.0,        0.0,        0.0, 2.0*L*L*Iz,
        -L*L*A,      0.0,       0.0,        0.0,        0.0,        0.0,  L*L*A,       0.0,      0.0,        0.0,        0.0,        0.0,
        0.0   , -12.0*Iz,       0.0,        0.0,        0.0,  -6.0*L*Iz,    0.0,   12.0*Iz,      0.0,        0.0,        0.0,  -6.0*L*Iz,
        0.0   ,      0.0,  -12.0*Iy,        0.0,   6.0*L*Iy,        0.0,    0.0,       0.0,  12.0*Iy,        0.0,   6.0*L*Iy,        0.0,
        0.0   ,      0.0,       0.0, -G*L*L*J/E,        0.0,        0.0,    0.0,       0.0,      0.0,  G*L*L*J/E,        0.0,        0.0,
        0.0   ,      0.0, -6.0*L*Iy,        0.0, 2.0*L*L*Iy,        0.0,    0.0,       0.0, 6.0*L*Iy,        0.0, 4.0*L*L*Iy,        0.0,
        0.0   , 6.0*L*Iz,       0.0,        0.0,        0.0, 2.0*L*L*Iz,    0.0, -6.0*L*Iz,      0.0,        0.0,        0.0, 4.0*L*L*Iz
    )*(E/L^3)
end

"""
Computes the mass matrix of the beam element
- A: Cross section.
- L: Length of the beam.
- ρ: Beam density.

Taken from:
Dynamic Analysis of Structures by John T. Katsikadelis: eq (11.9.11)
"""
function Mass_Matrix(A::Float64, L::Float64, ρ::Float64)::SMatrix{12, 12, Float64, 144}
    Ix::Float64 = 2.0*L*(A/π)*sqrt(A/π)/3.0 # 2*L*r^3/3 -> Second moment of area.
    r2::Float64 = Ix/A
    return SMatrix{12, 12, Float64, 144}(
        140.0,     0.0,     0.0,      0.0,      0.0,      0.0,  70.0,     0.0,     0.0,      0.0,      0.0,      0.0,
        0.0  ,   156.0,     0.0,      0.0,      0.0,   22.0*L,   0.0,    54.0,     0.0,      0.0,      0.0,  -13.0*L,
        0.0  ,     0.0,   156.0,      0.0,  -22.0*L,      0.0,   0.0,     0.0,    54.0,      0.0,   13.0*L,      0.0,
        0.0  ,     0.0,     0.0, 140.0*r2,      0.0,      0.0,   0.0,     0.0,     0.0,  70.0*r2,      0.0,      0.0,
        0.0  ,     0.0, -22.0*L,      0.0,  4.0*L*L,      0.0,   0.0,     0.0, -13.0*L,      0.0, -3.0*L*L,      0.0,
        0.0  ,  22.0*L,     0.0,      0.0,      0.0,  4.0*L*L,   0.0,  13.0*L,     0.0,      0.0,      0.0,  -3.0*L*L,
        70.0 ,     0.0,     0.0,      0.0,      0.0,      0.0, 140.0,     0.0,     0.0,      0.0,      0.0,      0.0,
        0.0  ,    54.0,     0.0,      0.0,      0.0,   13.0*L,   0.0,   156.0,     0.0,      0.0,      0.0,  -22.0*L,
        0.0  ,     0.0,    54.0,      0.0,  -13.0*L,      0.0,   0.0,     0.0,   156.0,      0.0,   22.0*L,      0.0,
        0.0  ,     0.0,     0.0,  70.0*r2,      0.0,      0.0,   0.0,     0.0,     0.0, 140.0*r2,      0.0,      0.0,
        0.0  ,     0.0,  13.0*L,      0.0, -3.0*L*L,      0.0,   0.0,     0.0,  22.0*L,      0.0,  4.0*L*L,      0.0,
        0.0  , -13.0*L,     0.0,      0.0,      0.0, -3.0*L*L,   0.0, -22.0*L,     0.0,      0.0,      0.0,  4.0*L*L,
    )*(ρ*A*L/420.0)
end