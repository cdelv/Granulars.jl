include("../../src/Granulars.jl")
l = 1.0
g = 1.0

function Calculate_Forces(particles::StructVector{Particle}, 
    neighborlist::Vector{Tuple{Int64, Int64, Float64}}, 
    conf::Config,
    kundall_particles::ExtendableSparseMatrix{Float64, Int64},
    kundall_walls::ExtendableSparseMatrix{Float64, Int64})
    q = particles.q[1]
    
    particles.τ[1] = particles.m[1]*g*l*SVector(
            2.0*(q.q2*q.q3 + q.q0*q.q1),
            -2.0*(q.q1*q.q3 - q.q0*q.q2), 
            0.0)
    
    return nothing
end

function Save_step(particles::StructVector{Particle}, file::String, i::Int, t::Float64)
    
    # Define The Header Depending On The Dimmension Of The Particles 
    H = ["t","x","y","z","vx","vy","vz","ax","ay","az","a1","a2","a3","wx","wy","wz","αx","αy","αz","m","Ix","Iy","Iz","rad"]

    # This syntax is a bit involved. We have a StrucArray.
    # For CSV.write to work, we need a vector of vectors and then create a table.
    # 1. Reinterpret gets rid of the SVectors on each row.
    # 2. Reshape creates the vector of vectors.
    # 3. transpose to print in the correct order.
    t = transpose(reshape(reinterpret(Float64,t*ones(length(particles))),(1,length(particles))))
    r = transpose(reshape(reinterpret(Float64,particles.r),(3,length(particles))))
    v = transpose(reshape(reinterpret(Float64,particles.v),(3,length(particles))))
    a = transpose(reshape(reinterpret(Float64,particles.a),(3,length(particles))))
    # q
    w = transpose(reshape(reinterpret(Float64,particles.w),(3,length(particles))))
    τ = transpose(reshape(reinterpret(Float64,particles.τ),(3,length(particles))))
    m = transpose(reshape(reinterpret(Float64,particles.m),(1,length(particles))))
    II = transpose(reshape(reinterpret(Float64,particles.I),(3,length(particles))))
    rad = transpose(reshape(reinterpret(Float64,particles.rad),(1,length(particles))))
    
    q = StructArray(quat_to_angle.(particles.q, :XYZ))
    a1 = transpose(reshape(reinterpret(Float64,q.a1),(1,length(particles))))
    a2 = transpose(reshape(reinterpret(Float64,q.a2),(1,length(particles))))
    a3 = transpose(reshape(reinterpret(Float64,q.a3),(1,length(particles))))

    # 4. The table creates the table to print.
    # Does allocations, how can I remove them?
    Table = CSV.Tables.table(hcat(t,r,v,a,a1,a2,a3,w,τ,m,II,rad))
    
    # 5. CSV.write allows the user to choose the buffer size. It does not have 
    # performance impact. The default is 4 Mb (2^22). We don't need to use that much memory.
    # CSV.write is the fastest writter that I could find.
    # bufsize=2^9 min for 3D.
    CSV.write(file*"_"*string(i)*".csv",Table, header=H, bufsize=2^9) 
    
    nothing
end

function main()
    W1 = Wall([1,0,0],[0,0,0])
    walls = [W1]
    p = Particle(w=[0, 1.0, 9.0])
    p = Set_I(p,SVector{3,Float64}([0.2, 0.2, 1.0]))
    p = Set_q(p,angle_to_quat(π/7, 0, -π/2,:XYZ))
    data = [p]
    conf = Config(0.5,0.0001,g=ones(3),walls=walls)
    Propagate(data, conf; vis_steps=20, file="Paraview/data", save=true)
end

main()