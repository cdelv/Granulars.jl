
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

function Check_Simulation(particles::StructVector{<:AbstractParticle})
    for i in eachindex(particles)
        @assert all(>(0.0), particles.I[i])
        @assert particles.m[i] > 0.0
        @assert particles.rad[i] > 0.0
        if !(norm(particles.q[i]) ≈ 1)
            particles.q[i] = particles.q[i]/norm(particles.q[i])
            println("Normalizing quaternion of particle ", i)
        end
    end
end

"""
Computes the inertia tensor of a spherical particle using the montecarlo integration.
For a sphere this is overkill but for the multispheres will be usefull. 

PUT THE EXACT INERTIA WITH WHEN THE MULTISPHERE IS DONE

https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/dd277ec654440f4c2b5b07d6c286c3fd_MIT16_07F09_Lec26.pdf
"""
function Compute_Inertia_Tensor(p::Particle, num::Int64=100000)::Matrix{Float64}
    Ixx::Float64 = 0.0
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
    return (p.m/T)*SMatrix{3,3}(Ixx, -Ixy, -Ixz, -Ixy, Iyy, -Ixy, -Ixy, Iyz, Izz)
end

"""
Returns the smallest angle between two vectors in radians where  0 <= angle(v, w) <= π
It uses atan(norm(cross(u,v)),dot(u,v)) insted of acos(dot(v,w) / (norm(v)*norm(w))) as it is 
more acurate. acos fails with small angles. 

Inspired by AngleBetweenVectors.jl.

https://people.eecs.berkeley.edu/~wkahan/MathH110/Cross.pdf (page 15)
https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab#answer_185622
"""
function angle(v::SVector{3},w::SVector{3})::Float64
    a::Float64 = atan(norm(cross(v,w)),dot(v,w))
end

# reference frame rotations
function Lab_to_body(v::SVector{3}, q::Quaternion{Float64})::SVector{3}
    return quat_to_dcm(q)*v
end

function Body_to_lab(v::SVector{3}, q::Quaternion{Float64})::SVector{3}
    return inv_rotation(quat_to_dcm(q))*v
end