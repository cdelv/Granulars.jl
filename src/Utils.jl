
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