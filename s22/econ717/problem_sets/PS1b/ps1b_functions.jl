#==
    This file defines all user-written functions called
    by edgel_ps1b.jl
==#


# load dependencies
using Parameters, Random, DataFrames, Distributions, LinearAlgebra, StatsBase

## Define the structure that contains "true" model parameters
@with_kw mutable struct Primitives

    # Unobeserved model primitives
    π₁::Float64 = 1.0  # wage in occupation 1
    π₂::Float64 = 1.0  # wage in occupation 2
    μ₁::Float64        # mean of ln(S₁) 
    μ₂::Float64        # mean of ln(S₂)
    σ₁::Float64        # std of ln(S₁)
    σ₂::Float64        # std of ln(S₂)
    ρ::Float64         # correlation between ln(S₁) and ln(S₂)

    # Secondary primitives to be filled in by Initialize()
    Σ::Array{Float64, 2}

end # Primitives struct 

### Define a function that initializes the primitives
function Initialize(; θ::Vector{Float64} = [1, 1, 0.5, 0.5, 0.25, 0.25, 0.4])

    # Unpack the parameters of the parameter vector
    π₁, π₂, μ₁, μ₂, σ₁, σ₂, ρ = θ

    # Initialize the primitives struct
    prim = Primitives(π₁=π₁, π₂=π₂, μ₁=μ₁, μ₂=μ₂, σ₁=σ₁, σ₂=σ₂,
        ρ=ρ, Σ=[σ₁^2 ρ*σ₁*σ₂; ρ*σ₁*σ₂ σ₂^2])

    # Return the primitives struct, parameter vector, and data
    return prim, θ

end # Initialize()

### Define a function that generates data, given primitives
function SimulateData(θ::Vector{Float64}, N::Int64; seed::Int64=115,
    w̲::Float64=-Inf)

    # set seed for reproducibility
    Random.seed!(seed)

    # unpack relevant primitives from θ
    π₁, π₂, μ₁, μ₂, σ₁, σ₂, ρ = θ
    Σ = [σ₁^2 ρ*σ₁*σ₂; ρ*σ₁*σ₂ σ₂^2]

    # generate an array of N simulated agents
    ε = rand(Distributions.MvNormal([0, 0], Σ), N)'
    w₁, w₂ = max.(π₁.*(μ₁ .+ ε[:, 1]), w̲), π₂.*(μ₂ .+ ε[:, 2])

    # return DataFrame of observed results 
    return DataFrame(i=1:N, j=1 .+ 1*(w₁ .< w₂), 
                        W=(w₁ .>= w₂).*w₁ + (w₁ .< w₂).*w₂)

end # SimulateData()

### Define the objective function for finding a parameter
### vector that results in p% of agents choosing 1
function OccChoiceObj(θ::Vector{Float64}; p::Float64 = 0.6, 
    N::Int64 = 1000, w̲::Float64 = -Inf, seed::Int64 = 3609)

    # return the squard difference between p and p̂
    return (p - mean(SimulateData(θ, N; w̲ = w̲, seed = seed).j .== 1))^2

end # OccChoiceObj()

### Define the SMM objective function for finding a parameter
function SMMObjFun(θ::Vector{Float64}, data::DataFrame; 
    N::Int64 = 1000, m::Int64 = length(θ), seed::Int64 = 115,
    W::Matrix{Float64} = 1.0*Matrix(I, m, m)) 

    # impose boundary conditions on k 
    @assert (m > 0) "m must be a positive integer"

    # simulate data
    sim = SimulateData(θ, N; seed = seed)

    # calculate simulated moments
    gθ = [
        mean(sim.j .== 1),
        mean(sim.W[sim.j .== 1]),
        mean(sim.W[sim.j .== 2]),
        var(sim.W[sim.j .== 1]),
        var(sim.W[sim.j .== 2]),
        skewness(sim.W[sim.j .== 1]),
        skewness(sim.W[sim.j .== 2])
    ]

    # calculate observed moments
    ĝ = [
        mean(data.j .== 1),
        mean(data.W[data.j .== 1]),
        mean(data.W[data.j .== 2]),
        var(data.W[data.j .== 1]),
        var(data.W[data.j .== 2]),
        skewness(data.W[data.j .== 1]),
        skewness(data.W[data.j .== 2])
    ]

    # subset to the first m moments
    gθ, ĝ = gθ[1:m], ĝ[1:m]

    # return the objective function
    return (gθ - ĝ)'*W*(gθ - ĝ)
end