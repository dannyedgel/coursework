#==
    This file defines the functions called by edgel_ps2.jl
==#

# Declare packages
using Distributions, Distributed, SharedArrays, StatsBase

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
    S = Array{Float64}(n, 2, s)
    S[:, 2, :] = rand(Bernoulli(param.p), (n, s))                            # X
    S[:, 1, :] = param.β + param.θ.*sim[:, 2, :] + rand(Normal(), (n, s))    # Y
    sim   = Simulation(S, s, n)

    return param, sim
end # initialize

# Function that computes θₙ given Y and X
function OLS(Y, X; modify = false)
    if modify
        return maximum(0.01, sum((X .- mean(X)).*Y)/sum((X .- mean(X))^2))
    else
        return sum((X .- mean(X)).*Y)/sum((X .- mean(X))^2)
    end
end # OLS

# Function that computes heteroskedasticity-robust asymptotic SE 
function AsympSE(Y, X, param::DGP)
    uhat = Y - param.β - param.θ.*X
    return sqrt((sum((X .- mean(X)).^2)^(-2))*sum(((X .- mean(X)).^2).*(uhat.^2)))
end # Asymptotic SE 

# Function that computes bootstrap standard errors
function BootstrapSE(Y, X; B::Int64 = 1000)

    # generate bootstrap samples of Y and X
    Ys, Xs = Y[rand(1:length(Y), (length(Y), B))], X[rand(1:length(Y), (length(Y), B))]

    # initialize vector of θhats
    @everywhere Θhat = SharedArray(B)

    # generate vector of θhats
    @async @distributed for b = 1:B
        Θhat[b] = OLS(Ys[:, b], X[:, b])
    end

    # return SE and IQR SE
    return (B - 1)^{-1}*sum((θhat .- mean(θhat).^2)), iqr(θhat) / (quantile(Normal(), .75) - quantile(Normal(), .25))
    
end # Bootstrap SE

# Function that computes score bootstrap standard errors
function ScoreBootstrap(Y, X; B::Int64 = 1000, Binomial = false)
    #== 
    TODO: use notes on score bootstrap to complete below
    # generate bootstrap samples of Y and X
    Ys, Xs = Y[rand(1:length(Y), (length(Y), B))], X[rand(1:length(Y), (length(Y), B))]

    # initialize vector of θhats
    @everywhere Θhat = SharedArray(B)

    # generate vector of θhats
    @async @distributed for b = 1:B
        Θhat[b] = OLS(Ys[:, b], X[:, b])
    end

    # return SE and IQR SE
    return (B - 1)^{-1}*sum((θhat .- mean(θhat).^2)), iqr(θhat) / (quantile(Normal(), .75) - quantile(Normal(), .25))
    ==#
end # Score Bootstrap

# Function that computes convergence rate
function ConvergenceRate(θhat, θ, se)

end # Convergence rate