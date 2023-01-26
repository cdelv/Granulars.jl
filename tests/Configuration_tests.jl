using Test
include("../src/Granulars.jl")

# CONSTRUCTORS
Wall(ones(SVector{3}),ones(SVector{3}))
W = Wall([1,1,1],[1,1,1])

Config(10.0,0.1,1.0e6,500.0,ones(SVector{3}),[W],900.0,0.4)
Config(10.0,0.1, walls=[W])
Config(10.0,0.1, walls=[W], g=[0,-9,0])