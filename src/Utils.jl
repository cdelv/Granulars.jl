
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