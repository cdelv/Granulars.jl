using StaticArrays
using DelimitedFiles
using LinearAlgebra
using Plots

const G=4*pi^2      #AU^3 Solar_mass^-1 year^-2

function main()
    dt=0.001      #0.001 years -> 9 hours 
    tf=100        #100 years
    vis_steps = 1000  #optional variable
    file = "data.csv" #optional variable
    bodies = start_from_file("S0stars.dat") #Initialize system
    convert_units(bodies, 8000, 1, 1) #Units on the file where on arcseconds, years and solar masses
    Propagate(bodies,tf,dt,vis_steps,file) 

    #Plot particles trayectory
    data=readdlm(file, ',')
    N = length(bodies)
    P = plot(data[:,2], data[:,3], data[:,4], label="SgrA*", zlabel="z [AU]")
    for i in 2:N
        plot!(data[:,2+7*(i-1)],data[:,3+7*(i-1)], data[:,4+7*(i-1)],label="S"*string(i-1))
    end
    xlabel!("x [AU]")
    ylabel!("y [AU]")
    savefig(P,"sagitarius_A*.png")

    #Animate particles trayectory
    trail=30
    anim = @animate for j in 1:size(data)[1]
        if j<trail+1
            e=1
        else
            e=j-trail
        end
        plot(data[e:j,2],data[e:j,3], data[e:j,4], lims = (-7e3, 8e3), label = "")
        scatter3d!((data[j,2],data[j,3], data[j,4]), label="SgrA*")
        for i in 2:N
            scatter3d!((data[j,2+7*(i-1)],data[j,3+7*(i-1)], data[j,4+7*(i-1)]), label="S"*string(i-1))
            plot!(data[e:j,2+7*(i-1)],data[e:j,3+7*(i-1)], data[e:j,4+7*(i-1)], label = "")
        end
    end
    gif(anim, "sagitarius_A*.gif")

    #random particles
    #=N=3
    particles = particle[]
    for i = 1:N
        rr = SVector(rand(), rand(), rand())
        vv =  SVector(rand(), rand(), rand())
        mass = rand()
        push!(particles, particle(rr, vv, mass))
    end=#
    
end

#Julia doesnt have classes, data structures instead
mutable struct particle
    r::SVector{3,Float64} #position
    v::SVector{3,Float64} #velocity
    a::SVector{3,Float64} #acceleration
    m::Float64            #masa
end
#new constructors
function particle(r,v,m) 
    return particle(r,v,SVector(0, 0, 0),m)
end
function particle(r,m)
    return particle(r,SVector(0, 0, 0),SVector(0, 0, 0),m)
end 

#particle methods
function clear_f(body)
    body.a=SVector(0, 0, 0)
end
function move_r(body,dt,cte=1)
    body.r+=body.v*dt*cte
end
function move_v(body,dt,cte=1)
    body.v+=body.a*dt*cte
end
function calc_f(body)
    clear_f.(body)
    for i in 1:length(body)
        for j in i+1:length(body)
            dr = body[j].r-body[i].r
            d = norm(dr)
            F = -G*dr/d^3
            body[i].a+=-F*body[j].m #3 Newton Law
            body[j].a+=F*body[i].m 
        end
    end
end

#start particles from data in file
function start_from_file(file)
    data=readdlm(file, comments=true)
    r = Array{Float64}[]
    v = Array{Float64}[]
    m = []
    for i in eachrow(data)
        push!(r,[i[1],i[2],i[3]])
        push!(v,[i[4],i[5],i[6]])
        append!(m,[i[7]])
    end
    return particle.(r,v,m) # . is the broadcast operation
end

#Evolve the system
function Propagate(data,tf,dt,vis_steps=10000,file="data.csv")
    N=trunc(Int, tf/dt) #number of steps
    t=0
    Save_step(t,data,file,true) #save initial condition

    for i in 1:N
        step_PEFRL(data,dt) #time integration
        t+=dt
        
        if i%vis_steps==0 #save data every vis_steps
            Save_step(t,data,file)
        end
    end
end
function Save_step(t,data,file,start=false)
    y = [t]
    for i in data
        y=vcat(y,i.r,i.v,i.m) #concatenate data in 1 array
    end
    
    if start
        f = open(file, "w")   #create or overwritte the file
        write(f, join(y,",")*"\n")
        close(f)
    else
        f = open(file, "a")  #append to the file
        write(f, join(y,",")*"\n")
        close(f)
    end
end
function step_PEFRL(body::Array{particle,1},Δt::Float64)
    const1 = 0.1786178958448091      #ζ
    const3 = -0.6626458266981849e-1  #χ
    const4 = -0.2123418310626054     #λ
    const2 = (1-2*const4)/2          #(1-2λ)/2
    const5 = 1-2*(const3+const1)     #1-2*(ζ+χ)
    
    for particle in body
        move_r(particle,Δt,const1)
        clear_f(particle)
    end
    calc_f(body)
    for particle in body
        move_v(particle,Δt,const2)
        move_r(particle,Δt,const3)
        clear_f(particle)
    end
    calc_f(body)
    for particle in body
        move_v(particle,Δt,const4)
        move_r(particle,Δt,const5)
        clear_f(particle)
    end
    calc_f(body)
    for particle in body
        move_v(particle,Δt,const4)
        move_r(particle,Δt,const3)
        clear_f(particle)
    end
    calc_f(body)
    for particle in body
        move_v(particle,Δt,const2)
        move_r(particle,Δt,const1)
    end
end
function convert_units(data, Rconv, Tconv, Mconv)
    for i in eachindex(data)
        data[i].r*=Rconv
        data[i].v*=Rconv/Tconv
        data[i].m*=Mconv
    end
end

main()