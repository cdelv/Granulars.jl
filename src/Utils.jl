"""
Damping coeficient as a function of the coeficient of restitution. Taken from:

An investigation of the comparative behaviour of alternative contact force models during inelastic collisions by Colin Thornton 

- en: Coeficient of restitution.
"""
function γ(en::Float64)::Float64
    if en >= 1.0
        return 0.0 
    end
    h1::Float64 = -6.918798
    h2::Float64 = -16.41105
    h3::Float64 = 146.8049
    h4::Float64 = -796.4559
    h5::Float64 = 2928.711
    h6::Float64 = -7206.864
    h7::Float64 = 11494.29
    h8::Float64 = -11342.18
    h9::Float64 = 6276.757
    h10::Float64 = -1489.915
    α::Float64 = en*(h1+en*(h2+en*(h3+en*(h4+en*(h5+en*(h6+en*(h7+en*(h8+en*(h9+en*h10)))))))))
    return sqrt((1.0/(1.0-(1.0+en)^2*exp(α)))-1.0)
end


"""
Returns a normalize vector unless the norm=0, then it returns [0,0,0].
- v: vector to normalize.
"""
function unitary(v::SVector{3, Float64})::SVector{3, Float64}
    if norm(v)==0
        return zeros(SVector{3})
    else 
        return normalize(v)
    end
end

"""
Returns a normalize quaternion. The norm has to be different than 0.
- q: quaternion to normalize.
"""
function unitary(q::Quaternion{Float64})::Quaternion{Float64}
    return q/norm(q)
end

"""
Computes the inertia tensor of a spherical particle using montecarlo integration.
For a sphere this is overkill but for the multispheres will be usefull. 
- p: particle to use in the calculations.
- num: number of samples for the integral. 

PUT THE EXACT INERTIA WHEN THE MULTISPHERE IS DONE

https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/dd277ec654440f4c2b5b07d6c286c3fd_MIT16_07F09_Lec26.pdf
"""
function Compute_Inertia_Tensor(p::Particle, num::Int64=50000)::Matrix{Float64}
    #=Ixx::Float64 = 0.0
    Iyy::Float64 = 0.0
    Izz::Float64 = 0.0
    Ixy::Float64 = 0.0
    Ixz::Float64 = 0.0
    Iyz::Float64 = 0.0
    T::Float64 = 0.0
    for row in eachrow(rand(Uniform(-p.rad,p.rad), num,3))
        if norm(row) <= p.rad
            @inbounds Ixx += row[2]^2+row[3]^2
            @inbounds Iyy += row[1]^2+row[3]^2
            @inbounds Izz += row[1]^2+row[2]^2
            @inbounds Ixy += row[1]*row[2]
            @inbounds Ixz += row[1]*row[3]
            @inbounds Iyz += row[2]*row[3]
            @inbounds T+=1.0
        end
    end
    return (p.m/T)*SMatrix{3,3}(Ixx, -Ixy, -Ixz, -Ixy, Iyy, -Ixy, -Ixy, Iyz, Izz)=#
    (2.0*p.m*p.rad*p.rad/5.0)*SMatrix{3,3}(1.0I(3))
end


"""
Sets the inertia tensor and orientation of a particle.
- p: particle to edit.
- num: number of samples for the integral. 

DOES ALLOCATIONS!

CHECK EIGENVECS
"""
function Set_Inertia(p::Particle, num::Int64=50000)::Particle
    Inertia::Matrix{Float64} = Compute_Inertia_Tensor(p::Particle, num)
    II::SVector{3, Float64} = SVector{3, Float64}(abs.(eigvals(Inertia))) # Principal axis inertia tensor (diagonal M).
    index::SVector{3, Int64} = sortperm(II, rev=true) # convention: biggest value in x, smallest in z.
    Axis::SMatrix{3,3, Float64} = SMatrix{3,3, Float64}(abs.(eigvecs(Inertia))[index,:]) # CHECK THIS!
    p = Set_I(p,II[index]) # sets the inertia in the principal axys
    p = Set_q(p,dcm_to_quat(DCM(transpose(Axis)))) # set the orientation of principal axis
    Set_q(p,p.q/norm(p.q)) # normalize the quaternion
end


"""
Returns the smallest angle between two vectors in radians where  0 <= angle(v, w) <= π
It uses atan(norm(cross(u,v)),dot(u,v)) insted of acos(dot(v,w) / (norm(v)*norm(w))) as it is 
more acurate. acos fails with small angles. 

Inspired by AngleBetweenVectors.jl.

https://people.eecs.berkeley.edu/~wkahan/MathH110/Cross.pdf (page 15)
https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab#answer_185622
"""
function angle(v::SVector{3, Float64},w::SVector{3, Float64})::Float64
    atan(norm(cross(v,w)),dot(v,w))
end


"""
Transforms a vector from the lab frame to the body frame.
- v: vector to transform.
- q: quaternion with the orientation of the body frame.
"""
function Lab_to_body(v::SVector{3, Float64}, q::Quaternion{Float64})::SVector{3, Float64}
    vect(inv(q)*v*q) # an other way of doing it
    #quat_to_dcm(q)*v
end

"""
Transforms a vector from the body frame to the lab frame.
- v: vector to transform.
- q: quaternion with the orientation of the body frame.
"""
function Body_to_lab(v::SVector{3, Float64}, q::Quaternion{Float64})::SVector{3, Float64}
    vect(q*v*inv(q)) # an other way of doing it
    #inv_rotation(quat_to_dcm(q))*v
end


"""
Compute the time-derivative of the quaternion `qba` that rotates a reference
frame `a` into alignment to the reference frame `b` in which the angular
velocity of `b` with respect to `a`, and represented in `b`, is `wba_b`.

This is a re rwitte from ReferenceFrameRotations.jl wich makes the computation faster.
The original function has some stuff I dont need. 
"""
function dquat(qba::Quaternion{Float64}, w::SVector{3, Float64})::Quaternion{Float64}
    # Return the time-derivative.
    #         1         x
    #   dq = --- (wba_b) . q
    #         2
    @inbounds Quaternion(
                          - w[1] / 2 * qba.q1 - w[2] / 2 * qba.q2 - w[3] / 2 * qba.q3,
        w[1] / 2 * qba.q0                     + w[3] / 2 * qba.q2 - w[2] / 2 * qba.q3,
        w[2] / 2 * qba.q0 - w[3] / 2 * qba.q1                     + w[1] / 2 * qba.q3,
        w[3] / 2 * qba.q0 + w[2] / 2 * qba.q1 - w[1] / 2 * qba.q2
    )
end

"""
Calculates the beam orientation. The beam director vector must be aligned with the x axis 
in the beam reference frame. 
- n: Beam director vector. 
"""
function Beam_Orientation(n::SVector{3, Float64})::Quaternion{Float64}
    ex::SVector{3} = SVector(1.0,0.0,0.0) # x axis
    Angle::Float64 = angle(ex,n) # rotation angle
    u::SVector{3} = cross(ex,n)/(sin(Angle)+1e-12) # rotation axis
    unitary(angleaxis_to_quat(EulerAngleAxis(Angle, u)))
end

"""
Estimates the correct time step for the simulation using a percentage of the PWave timestep (usually 50-80% of it)
https://gitlab.com/yade-dev/trunk/-/blob/master/preprocessing/dem/Shop_01.cpp#L469-486

"""
function PWaveTimeStep(particles::Vector{<:Particle})::Float64
    dt::Float64 = Inf

    for i in particles
        dt = min(dt, i.rad*sqrt(Get_Density(i)/i.E))
    end

    return dt
end


"""
Creates the walls needed for a box with a corner in (0,0,0).
- lx: Dimension of the box in the x direction.
- ly: Dimension of the box in the y direction.
- lz: Dimension of the box in the z direction.
- E: Young modulus for the walls.
- G: Shear modulus for the walls.
"""
function Create_Box(lx::Real, ly::Real, lz::Real; E::Real=1.0e6, ν::Real=0.2)::Vector{Wall}
    Lx::Float64 = Float64(lx)
    Ly::Float64 = Float64(ly)
    Lz::Float64 = Float64(lz)

    # X coordinate walls
    W1::Wall = Wall(SVector(1.0,0.0,0.0), SVector(0.0,0.0,0.0), E=E, ν=ν)
    W2::Wall = Wall(SVector(-1.0,0.0,0.0),SVector(Lx,0.0,0.0), E=E, ν=ν)
    
    # Y coordinate walls
    W3::Wall = Wall(SVector(0.0,1.0,0.0), SVector(0.0,0.0,0.0), E=E, ν=ν)
    W4::Wall = Wall(SVector(0.0,-1.0,0.0),SVector(0.0,Ly,0.0), E=E, ν=ν)
    
    # Z coordinate walls
    W5::Wall = Wall(SVector(0.0,0.0,1.0), SVector(0.0,0.0,0.0), E=E, ν=ν)
    W6::Wall = Wall(SVector(0.0,0.0,-1.0),SVector(0.0,0.0,Lz), E=E, ν=ν)
    
    [W1,W2,W3,W4,W5,W6]
end


"""

"""
function Actions_Before_Time_Step!(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    friction_spring_particles::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    friction_spring_walls::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    beam_bonds::Dict{Tuple{Int64, Int64}, Int64},
    beams::StructVector{Beam},
    fixed_spheres::Vector{Int64},
    static::Bool,
    t::Float64)
    nothing
end

"""

"""
function Actions_After_Time_Step!(particles::StructVector{<:AbstractParticle}, 
    conf::Config, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}},
    friction_spring_particles::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    friction_spring_walls::Dict{Tuple{Int64, Int64}, SVector{3, Float64}},
    beam_bonds::Dict{Tuple{Int64, Int64}, Int64},
    beams::StructVector{Beam},
    fixed_spheres::Vector{Int64},
    static::Bool,
    t::Float64)
    nothing
end

"""
"""
function sparcify(v::Float64, eps::Float64=5.0e-14)::Float64
    if abs(v) < eps
        return 0.0
    else
        return v
    end
end

function CalculateA(target_e::Float64, target_v::Float64, Eeff::Float64, Reff::Float64, meff::Float64)::Float64

    function Target_Epsilon(v_star::Float64)::Float64
        Epsilon(v_star)-target_e
    end

    target_v_star::Float64 = find_zero((Target_Epsilon, dEpsilon_dv), 1.0, Roots.Newton())

    get_A(target_v, target_v_star, Eeff, Reff, meff)
end

function Epsilon(v_star::Float64)::Float64
    # 3/6 coefficients
    a_i::SVector{4, Float64} = SVector(1.0, 1.07232, 0.574198, 0.141552)
    b_i::SVector{7, Float64} = SVector(1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585)
    
    # Initialize sum to 0
    A::Float64 = 0.0
    B::Float64 = 0.0
    n::Float64 = 0.0

    for i in a_i
        A += i*v_star^n
        n += 1.0
    end
    
    n = 0.0
    for i in b_i
        B += i*v_star^n
        n += 1.0
    end
    
    return A/B
end

function dEpsilon_dv(v_star::Float64)::Float64
    # 3/6 coefficients
    a_i::SVector{4, Float64} = SVector(1.0, 1.07232, 0.574198, 0.141552)
    b_i::SVector{7, Float64} = SVector(1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585)
    
    # Initialize sum to 0
    A::Float64 = 0.0
    B::Float64 = 0.0
    dA::Float64 = 0.0
    dB::Float64 = 0.0
    n::Float64 = 0.0

    for i in a_i
        A += i*v_star^n
        n += 1.0
    end
    
    n = 0.0
    for i in b_i
        B += i*v_star^n
        n += 1.0
    end
    
    # Derivative
    da_i::SVector{3, Float64} = SVector(1.07232, 0.574198, 0.141552)
    db_i::SVector{6, Float64} = SVector(1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585)
    
    n = 0.0
    for i in da_i
        dA += (n+1)*i*v_star^n
        n += 1.0
    end
    
    n = 0.0
    for i in db_i
        dB += (n+1)*i*v_star^n
        n += 1.0
    end
    
    return (dA*B - A*dB)/(B*B)
end

function get_A(v::Float64, v_star::Float64, Eeff::Float64, Reff::Float64, meff::Float64)::Float64
    rho::Float64 = 4.0*Eeff*sqrt(Reff)/3.0
    beta::Float64 = v_star*v_star*v^(-0.2)      # ^-1/5 = ^-2/10
    2.0*beta*(rho/meff)^(-0.4)/3.0              # ^-2/5
end
        
function get_v_star(v::Float64, A::Float64, Eeff::Float64, Reff::Float64, meff::Float64)::Float64
    rho::Float64 = 4.0*Eeff*sqrt(Reff)/3.0
    beta::Float64 = 3.0*A*(rho/meff)^(0.4)/2.0; # ^2/5
    sqrt(beta)*v^(0.1)                          # ^1/10
end 