"""
Prints the information of all particles to a .csv file.
CSV.write is the fastest writer I could find. 
- TO DO: PRINT WALLS INFORMATION, CHECK IF THE DIRECTORY EXISTS, AND DELETE OLD DATA. 
- particles: StructArray of particles.
- file: Name of the file group to save all the csvs.
- i: Integer that defines the number of frame to print. 
"""
function Save_step(particles::StructVector{Particle{D}}, file::String, i::Int) where D
    
	# Define The Header Depending On The Dimmension Of The Particles 
    if D==2
        H = ["x","y","vx","vy","ax","ay","m","rad"]
    elseif D==3
        H = ["x","y","z","vx","vy","vz","ax","ay","vz","m","rad"]
    end
    
    m::Int32 = length(H)
    
    # This syntax is a bit involved, but it's to avoid allocations. We have a StrucArray.
    # For CSV.write to work, we need a vector of vectors and then create a table.
    # 1. Reinterpret gets rid of the SVectors on each row
    # 2. Reshape creates the vector of vectors.
    # 3. transpose to print in the correct order.
    # 4. The table creates the table to print.
    # 5. CSV.write allows the user to choose the buffer size. It does not have 
    # performance impact. The default is 4 Mb (2^22), so we don't need to use that much memory.
    # Creating the table is very fast. 
    CSV.write(file*"_"*string(i)*".csv",CSV.Tables.table(transpose(reshape(reinterpret(Float64,particles),(m,length(particles))))), header=H, bufsize=2^9) #bufsize=2^9 min for 3D. 
    
    nothing
end