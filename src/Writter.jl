"""
Print the information of all particles to a .csv
CSV.write if the fastest writter I could find. 
TO DO: PRINT WALLS INFORMATION, CHECH IF DIRECTORY EXISTS, DELETE OLD DATA. 
FINISH DESCRIPTION
"""
function Save_step(particles::StructVector{Particle{D}}, file::String, i::Int) where D
    if Get_Dim(particles[1])==2
        H = ["x","y","vx","vy","ax","ay","m","rad"]
    elseif Get_Dim(particles[1]) == 3
        H = ["x","y","z","vx","vy","vz","ax","ay","vz","m","rad"]
    end
    
    m::Int32 = Get_Dim(particles[1])*3 + 2
    
    CSV.write(file*"_"*string(i)*".csv",CSV.Tables.table(transpose(reshape(reinterpret(Float64,particles),(m,length(particles))))), 
        header=H, bufsize=2^9) #bufsize=2^9 min for 3D no cambia el performance, solo la memoria alocada. 
    nothing
end