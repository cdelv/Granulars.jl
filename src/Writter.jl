"""
Prints the information of all particles to a .csv file.
CSV.write is the fastest writer I could find. 
- TO DO: PRINT WALLS INFORMATION, CHECK IF THE DIRECTORY EXISTS, AND DELETE OLD DATA. 

- particles: StructArray of particles.
- file: Name of the file group to save all the csvs.
- i: Integer that defines the number of the frame to print.
- t: Simulation time. 
- rot_seq: Rotation sequence to use for the angles.
"""
function Save_step(particles::StructVector{Particle}, file::String, i::Int, t::Float64, rot_seq::Symbol=:XYZ)
    # Define The Header of the CSV file
    H = ["ID","t","x","y","z","vx","vy","vz","ax","ay","az","a1","a2","a3","wx","wy","wz","τx","τy","τz","m","Ix","Iy","Iz","rad"]

    # Arrange data for printing
    data = Proces_Particles.(1:length(particles), t, particles, rot_seq)

    # Create the table to print.
    Table = CSV.Tables.table(reduce(hcat, data)')
    
    # 5. CSV.write allows the user to choose the buffer size. It does not have 
    # performance impact. The default is 4 Mb (2^22). We don't need to use that much memory.
    # CSV.write is the fastest writter that I could find.
    # bufsize=2^9 min for 3D.
    CSV.write(file*"_"*string(i)*".csv",Table, header=H, bufsize=2^9) 
    
    nothing
end


"""
Rearranges the information on a particle to make it suitable for printing.
- i: Index of the particle.
- t: Simulation time. 
- p: Particle to print.
- rot_seq: Rotation sequence to use for the angles.
"""
function Proces_Particles(i::Int64, t::Float64, p::Particle, rot_seq::Symbol)::Vector{Float64}
    ϕ::EulerAngles{Float64} = quat_to_angle(p.q, rot_seq)
    w::SVector{3, Float64} = Body_to_lab(p.w,p.q)
    τ::SVector{3, Float64} = Body_to_lab(p.w,p.q)
    @inbounds [
        i, t, p.r[1],p.r[2],p.r[3],p.v[1],p.v[2],p.v[3],p.a[1],p.a[2],p.a[3],
        ϕ.a1,ϕ.a2,ϕ.a3,w[1],w[2],w[3],τ[1],τ[2],τ[3],
        p.m, p.I[1], p.I[2], p.I[3], p.rad
        ]
end