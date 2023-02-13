
"""
Structure that stores information of the beams. The beams are assumed cilindrical. 
- E: Young Modulus.
- G: Shear Modulus. 
- J: Torsion constant.
- I: Second area moment. 
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

    # Angle displacement
    Δq::Quaternion{Float64} = unitary(qbj ∘ inv(qbi))
    Δϕ::EulerAngles{Float64} = quat_to_angle(Δq, :XYZ)      

    # Create the 12x12 transformation matrix and calculate forces and torques.
    #@inbounds Δs::SVector{12,Float64} = 0.5*SVector(Δr[1], Δr[2], Δr[3], Δϕ.a1, Δϕ.a2, Δϕ.a3, -Δr[1], -Δr[2], -Δr[3], -Δϕ.a1, -Δϕ.a2, -Δϕ.a3)
    #F::SVector{12,Float64} = Stifness_Matrix(beams.A[k], beams.L[k], beams.E[k], beams.G[k], beams.I[k], beams.J[k])*Δs

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

    # Damping force according to Raleigh damping model BROKEN!!!
    vi = Lab_to_body(particles.v[i], qbi)
    vj = Lab_to_body(particles.v[j], qbi)
    wi = Lab_to_body(Body_to_lab(particles.w[i],particles.q[i]), qbi)
    wj = Lab_to_body(Body_to_lab(particles.w[j],particles.q[j]), qbi)
    @inbounds V::SVector{12,Float64} = (vi[1], vi[2], vi[3], wi[1], wi[2], wi[3], vj[1], vj[2], vj[3], wj[1], wj[2], wj[3])
    F -= (0.0*Stifness_Matrix(beams.A[k], beams.L[k], beams.E[k], beams.G[k], beams.I[k], beams.J[k]) + 1.0*Mass_Matrix(beams.A[k], beams.L[k], particles.m[i], particles.rad[i]))*V

    # Add forces and torques to the particles.
    @inbounds particles.a[i] += Body_to_lab(SVector(F[1], F[2], F[3]), qbi)/particles.m[i] 
    @inbounds particles.τ[i] += Lab_to_body(Body_to_lab(SVector(F[4], F[5], F[6]), qbi), particles.q[i])

    @inbounds particles.a[j] += Body_to_lab(SVector(F[7], F[8], F[9]), qbi)/particles.m[j]
    @inbounds particles.τ[j] += Lab_to_body(Body_to_lab(SVector(F[10],F[11],F[12]), qbi), particles.q[j])
    
    return nothing
end


"""
Computes the stifness matrix of the beam element
- A: Cross section.
- L: Length of the beam.
- E: Young Modulus.
- G: Shear Modulus.
- J: Torsion constant.  
"""
function Stifness_Matrix(A::Float64, L::Float64, E::Float64, G::Float64, I::Float64, J::Float64)::SMatrix{12, 12, Float64, 144}
    # we assume a cilindrical beam
    Iy::Float64 = I      # π*r^4/4 -> Second area moment.
    Iz::Float64 = I       # π*r^4/4 -> Second area moment.
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
Computes the mass matrix of the beam element
- A: Cross section.
- L: Length of the beam.
- m: Mass of one of the particles (for computing the density)
- rad: Radius of one of the particles (for computing the density)

Taken from:
http://what-when-how.com/the-finite-element-method/fem-for-frames-finite-element-method-part-1/
"""
function Mass_Matrix(A::Float64, L::Float64, m::Float64, rad::Float64)::SMatrix{12, 12, Float64, 144}
    a::Float64 = L/2.0
    Ix::Float64 = 2.0*L*(A/π)*sqrt(A/π)/3.0 # 2*L*r^3/3 -> Second area moment.
    rx::Float64 = sqrt(Ix/A)
    ρ::Float64 = 3.0*m/(4.0*π*rad^3)
    return SMatrix{12,12}(
        70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 35.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 78.0, 0.0, 0.0, 0.0, 22.0*a, 0.0, 27.0, 0.0, 0.0, 0.0, -13.0*a,
        0.0, 0.0, 78.0, 0.0, -22.0*a, 0.0, 0.0, 0.0, 27.0, 0.0, 13.0*a, 0.0,
        0.0, 0.0, 0.0, 70.0*rx*rx, 0.0, 0.0, 0.0, 0.0, 0.0, -35.0*rx*rx, 0.0, 0.0,
        0.0, 0.0, -22.0*a, 0.0, 8.0*a*a, 0.0, 0.0, 0.0, -13.0*a, 0.0, -6.0*a*a, 0.0,
        0.0, 22*a, 0.0, 0.0, 0.0, 8.0*a*a, 0.0, 13.0*a, 0.0, 0.0, 0.0, -6.0*a*a,
        35.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 27.0, 0.0, 0.0, 0.0, 13.0*a, 0.0, 78.0, 0.0, 0.0, 0.0, -22.0*a,
        0.0, 0.0, 27.0, 0.0, -13.0*a, 0.0, 0.0, 0.0, 78.0, 0.0, 22.0*a, 0.0,
        0.0, 0.0, 0.0, -35.0*rx*rx, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0*rx*rx, 0.0, 0.0,
        0.0, 0.0, 13.0*a, 0.0, -6.0*a*a, 0.0, 0.0, 0.0, 22.0*a, 0.0, 8.0*a*a, 0.0,
        0.0, -13.0*a, 0.0, 0.0, 0.0, -6.0*a*a, 0.0, -22.0*a, 0.0, 0.0, 0.0, 8.0*a*a
    )*(ρ*A*a/105.0)
end