
"""
Structure that stores information of the beams
"""
struct Beam
    r0::SVector{3, Float64}
    θ0::SVector{3, Float64}

    E::Float64
    G::Float64

    r::Float64
    L::Float64
end


"""
Computes the stifness matrix of the beam element
- r: Radius of the cross section.
- L: Length of the beam.
- E: Young Modulus.
- G: Shear Modulus. 
"""
function T_beam(r::Float64,L::Float64,E::Float64,G::Float64)::SMatrix{12, 12, Float64, 144}
    # we assume cilindrical beam
    A::Float64 = 2*π*r*r # cross section.
    J::Float64 = π*r^4/2 # Torsion constant.
    Iy::Float64 = π*r^4/4 # Second area moment.
    Iz::Float64 = Iy # Second area moment.
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

# Compute beam for

# Create beams

# Add beams