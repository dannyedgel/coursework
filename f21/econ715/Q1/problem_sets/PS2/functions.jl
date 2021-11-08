#==
    This file defines the functions called by edgel_ps2.jl
==#

# Declare packages
using Distributions, StatsBase, LinearAlgebra, 
        Statistics, Parameters

# Define the parameters of the data generation process
@with_kw mutable struct DGP
    β       ::Float64   = 0.0
    θ       ::Float64   = 1.0
    p       ::Float64   = 0.1  # Bernouli parameter
end # data generating process

# Simulations structure
mutable struct Simulation 
    S::Array{Float64}
    s::Int64 
    n::Int64
end # simulation

# Function that initializes s simulations of n observations
# from the DGP, reuturning an n by 2 by s array
function Initialize(n::Int64; s::Int64 = 1)
    param = DGP()
    # generate simulated data
    S = zeros(n, 2, s)
    S[:, 2, :] = rand(Bernoulli(param.p), (n, s))                          # X
    S[:, 1, :] = param.β .+ param.θ.*S[:, 2, :] .+ rand(Normal(), (n, s))    # Y
    sim   = Simulation(S, s, n)

    return param, sim
end # initialize

# Function that computes θₙ given Y and X
function OLS(Y, X; modify = false, returnβ = false)
    X = [ones(length(X), 1) X] # add intercept to X
    β = (X'*X) \ (X'*Y)

    if modify
        β[2] = maximum([0.01 β[2]])
    end

    if returnβ
        return [β[1]; β[2]]
    else
        return β[2]
    end
end # OLS

# Function that computes heteroskedasticity-robust asymptotic SE 
function AsympSE(Y, X)
    uhat = Y - [ones(length(X), 1) X]*OLS(Y, X; returnβ = true)
    return sqrt((sum((X .- mean(X)).^2)^(-2))*sum(((X .- mean(X)).^2).*(uhat.^2)))
end # Asymptotic SE 

# Function that computes bootstrap standard errors
function BootstrapSE(Y, X; B::Int64 = 1000, modify = false)

    # generate bootstrap samples of Y and X
    Ys, Xs = Y[rand(1:length(Y), (length(Y), B))], X[rand(1:length(Y), (length(Y), B))]

    # initialize vector of θhats
    Θhat = zeros(B)

    # generate vector of θhats
    for b = 1:B
        Θhat[b] = OLS(Ys[:, b], Xs[:, b]; modify = modify)
    end
    Θhat = filter(!isnan, Θhat)
    B    = length(Θhat)

    # return SE and IQR SE
    return ((B-1)^(-1))*sum((Θhat .- mean(Θhat)).^2), iqr(Θhat) / (quantile(Normal(), .75) - quantile(Normal(), .25))
    
end # Bootstrap SE

# Function that computes score bootstrap standard errors
function ScoreBootstrap(Y, X; B::Int64 = 1000, binomial = false)

    # generate vector of random shocks
    if binomial
        ε = 2*rand(Binomial(), (length(Y), B))
    else
        ε = rand(Normal(1, 1), (length(Y), B))
    end

    # initialize vector of θcrosses
    θcross, β = zeros(B), OLS(Y, X; returnβ = true)
    uhat      = Y - [ones(length(X), 1) X]*β

    # generate vector of θcrosses
    for b = 1:B
        θcross[b] = β[2] + sum((X .- mean(X)).*uhat.*ε[:, b])/sum((X .- mean(X)).^2)
    end

    # return SE
    return ((B-1)^(-1))*sum((θcross .- mean(θcross)).^2)

    # sqrt((length(Y)^(-1))*(sum((X .- mean(X)).*uhat)/sum((X .- mean(X)).^2))^2)

end # Score Bootstrap

# Function that outputs mean, sd, and coverage rate
function seStats(θhat, se, param::DGP)
    return mean(filter(!isnan,se)), std(filter(!isnan,se)), sum((param.θ .>= (θhat .- 1.96.*se)) .& (param.θ .<= (θhat .+ 1.96.*se)))/length(filter(!isnan,se))
end # Convergence rate