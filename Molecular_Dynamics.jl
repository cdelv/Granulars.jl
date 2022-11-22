#Usefull macros
# @time @btime @benchmark @allocated @code_warntype @simd @turbo @inbounds

using StaticArrays, Plots, BenchmarkTools, LinearAlgebra, DelimitedFiles

#PEFRL constants
const const1 = 0.1644986515575760     #ζ
const const3 = 0.1235692651138917e1   #χ
const const4 = -0.2094333910398989e-1 #λ
const const2 = (1-2*const4)/2         #(1-2λ)/2
const const5 = 1-2*(const3+const1);   #1-2*(ζ+χ)

#Physical constants and parameters
struct Config
    g::Float64
    K::Float64
    dt::Float64
    Nwalls::Int
end

#Atom data Type, general for any dimension
struct Atom{D}
    r::SVector{D, Float64}
    v::SVector{D, Float64}
    a::SVector{D, Float64}
    m::Float64
    rad::Float64
end

#=
Constructors are not type-stable and slow
Allow for convenient Atom declarations
For speed, use the default one
=#
function Atom(r::Array{<:Real,1},v::Array{<:Real,1},a::Array{<:Real,1},m::Real,rad::Real)::Atom
    n = length(r)
    return Atom(SVector{n,Float64}(r),SVector{n,Float64}(v),SVector{n,Float64}(a),Float64(m),Float64(rad))
end
function Atom(r::Array{<:Real,1},v::Array{<:Real,1},m::Real,rad::Real)::Atom
    n = length(r)
    return Atom(SVector{n,Float64}(r),SVector{n,Float64}(v),SVector{n,Float64}(zeros(n)),Float64(m),Float64(rad))
end
function Atom(r::Array{<:Real,1},m::Real,rad::Real)::Atom
    n = length(r)
    return Atom(SVector{n,Float64}(r),SVector{n,Float64}(zeros(n)),SVector{n,Float64}(zeros(n)),Float64(m),Float64(rad))
end

#=
Atom Methods, all of them are type stable, 1 allocation 
O(1), fast ~ 20 ns
=#
function Move_r(A::Atom{T}, dt::Float64, cte=1.0::Float64)::Atom{T} where {T}
    return Atom(A.r+A.v*dt*cte,A.v,A.a,A.m,A.rad)
end
function Move_v(A::Atom{T}, dt::Float64, cte=1.0::Float64)::Atom{T} where {T}
    return Atom(A.r,A.v+A.a*dt*cte,A.a,A.m,A.rad)
end
function Add_a(A::Atom{T}, a::SVector)::Atom{T} where {T}
    return Atom(A.r,A.v,A.a+a,A.m,A.rad)
end
function Set_a(A::Atom{T}, a::SVector)::Atom{T} where {T}
    return Atom(A.r,A.v,a,A.m,A.rad)
end
function Clear_a(A::Atom{T})::Atom{T} where {T}
    return Atom(A.r,A.v,zeros(SVector{length(A.a)}),A.m,A.rad)
end

#=
Type stable, 0 allocations, O(N²)
=#
function Calculate_Forces(data::T, c::Config) where {T}
    n::Int = length(data)
    g = MArray(0*data[1].r); g[2]=-c.g #Gravity goes in the 2 coordinate direction
    
    #Clear forces and add gravity
    for i in 1:n
        @inbounds data[i]=Set_a(data[i],SVector(g))
    end
    
    #Calculate interactions
    for i in 1:n
        #Pair interactions
        for j in i+1:n
            rij = @inbounds data[j].r-data[i].r
            s = @inbounds data[i].rad+data[j].rad-norm(rij)
            
            #Check if there is an interaction
            if s<0.0; continue; end
            
            F = c.K*s^1.5*rij/norm(rij)
            @inbounds data[i]=Add_a(data[i],-F/data[i].m)
            @inbounds data[j]=Add_a(data[j], F/data[j].m)
        end
    end
end

#=
Type stable, 0 allocations, O(N)
Without Calculate_Forces ~ 20N ns, N size of the data
Supports different types of arrays
=#
function PEFRL_time_step(data::T, c::Config) where {T}
    n::Int = length(data)
    for i in c.Nwalls+1:n                      
        @inbounds data[i]=Move_r(data[i],c.dt,const1) 
    end
    Calculate_Forces(data,c)
    for i in c.Nwalls+1:n                     
        @inbounds data[i]=Move_v(data[i],c.dt,const2) 
        @inbounds data[i]=Move_r(data[i],c.dt,const3)
    end
    Calculate_Forces(data,c)
    for i in c.Nwalls+1:n
        @inbounds data[i]=Move_v(data[i],c.dt,const4)
        @inbounds data[i]=Move_r(data[i],c.dt,const5)
    end
    Calculate_Forces(data,c)
    for i in c.Nwalls+1:n
        @inbounds data[i]=Move_v(data[i],c.dt,const4)
        @inbounds data[i]=Move_r(data[i],c.dt,const3)
    end
    Calculate_Forces(data,c)
    for i in c.Nwalls+1:n
        @inbounds data[i]= Move_v(data[i],c.dt,const2)
        @inbounds data[i]= Move_r(data[i],c.dt,const1)
    end
end

#=
Type stable
=#
function Propagate(data::T, c::Config, tf::Float64; vis_steps=1000::Int, file="data.csv"::String) where {T}
    N::Int = trunc(Int, tf/c.dt) #number of steps
    t::Float64 = 0.0
    Save_step(t,data,file,true) #save initial condition

    for i in 1:N
        PEFRL_time_step(data,c) #temporal integration
        t+=c.dt
        
        if i%vis_steps==0 #save data every vis_steps
            Save_step(t,data,file)
        end
    end
end

#=
Not type stable but kinda fast
=#
function Save_step(t::Float64, data::T, file::String, start=false::Bool) where {T}
    y::String = string(t)

    for i in 1:length(data)
        y=y*","*join(data[i].r,",")*","
        y=y*join(data[i].v,",")*","
        y=y*join(data[i].m,",")*","
        y=y*join(data[i].rad,",")
    end
    
    if start
        f = open(file, "w")   #Create or overwrite the file
        write(f, y*"\n")
        close(f)
    else
        f = open(file, "a")   #Append to the file
        write(f, y*"\n")
        close(f)
    end
end