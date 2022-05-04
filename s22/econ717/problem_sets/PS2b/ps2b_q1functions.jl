#==
    This file defines all user-written functions called
    by edgel_ps2b.jl for question 1
==#


# load dependencies
using Parameters, Random, DataFrames, Distributions, LinearAlgebra, StatsBase

## Define the structure that contains "true" model parameters
@with_kw mutable struct Primitives

    # DGP for hourly earnings
    lnwᵢ::Function = (s, a, ϵ) -> 1 + 0.05*s + 0.1*a + ϵ
    ϵDist::Normal{Float64} = Distributions.Normal(0, 0.5)
    aDist::Normal{Float64} = Distributions.Normal(0, 4)

    # DGP for units of schooling
    sᵢ::Function = (a, z₁, z₂, η) -> 3*a + z₁ + z₂ + η
    ηDist::Normal{Float64}   = Distributions.Normal(0, 1)
    z₁Dist::Normal{Float64}  = Distributions.Normal(0, 0.1)
    z₂Dist::Normal{Float64}  = Distributions.Normal(0, 25)
    z₃Dist::Uniform{Float64} = Distributions.Uniform(0, 25)

end # Primitives struct 

### Define a function that generates data, given primitives
function SimulateData(N::Int64; seed::Int64=115)

    # unpack model primitives 
    prim = Primitives()
    @unpack aDist, ϵDist, sᵢ, ηDist, z₁Dist, z₂Dist, z₃Dist, lnwᵢ = prim

    # set seed for reproducibility
    Random.seed!(seed)

    # generate an array of N simulated agents
    a, η, ϵ = rand(aDist, N)', rand(ηDist, N)', rand(ϵDist, N)'
    z₁, z₂, z₃ = rand(z₁Dist, N)', rand(z₂Dist, N)', rand(z₃Dist, N)'


    # return matrix of observed results 
    return [lnwᵢ.(sᵢ(a, z₁, z₂, η), a, ϵ)' sᵢ.(a, z₁, z₂, η)' z₁' z₂' z₃']

end # SimulateData()

### Define an OLS function
function OLS(Y::Vector{Float64}, X::Array{Float64}; 
    intercept::Bool = true)

    # if intercept is specified, add a column of ones to X
    if intercept; X = [ones(size(X, 1), 1) X]; end

    # calculate coefficients 
    β = inv(X'X)*X'Y

    # calculate variance-covariance matrix
    Ω = (Y - X*β)'*(Y - X*β)/(size(Y, 1) - size(X, 2))
    V = inv(X'X).*Ω

    # return the point estimates and variance-covariance matrix
    return β, V

end # OLS()

### Define a 2SLS function 
function TwoSLS(Y::Vector{Float64}, X::Array{Float64},
    Z::Array{Float64}; intercept::Bool=true, Ftest::Bool=false,
    endogIndex::Int64=1)

    # if intercept is specified, add a column of ones to X
    if intercept
        X = [ones(size(X, 1), 1) X]
        endogIndex += 1
    end

    # adjust Z to include all exogenous variables
    Z = [Z X[:, setdiff(1:size(X, 2), endogIndex)]]

    # calculate coefficients
    β = inv(X'Z * inv(Z'Z) * Z'X) * X'Z * inv(Z'Z) * Z'Y

    # calculate variance-covariance matrix
    #n = size(Y, 1)
    #Qₓ, Qᵥ = (1 / n) * X'Z, (1 / n) * Z'Z
    #Ω = (1 / n) * sum(sum(Z .^ 2, dims=2) .* ((Y - X * β) .^ 2))
    #V = inv(Qₓ * inv(Qᵥ) * Qₓ) * (Qₓ * inv(Qᵥ) * Ω * inv(Qᵥ) * Qₓ) * inv(Qₓ * inv(Qᵥ) * Qₓ)
    Ω = (Y - X*β)'*(Y - X*β)/(size(Y, 1) - size(X, 2))
    V = inv(X'Z*inv(Z'Z)*Z'X).*Ω 

    # if specified, calculate the F statistic from the 
    # first-stage regression
    if Ftest
        Xₑ = X[:, endogIndex]
        Z₂ = X[:, setdiff(1:size(X, 2), endogIndex)]

        b, v = OLS(Xₑ, Z; intercept=false)
        σ¹ = (Xₑ - Z*b)' * (Xₑ - Z*b)

        b, v = OLS(Xₑ, Z₂; intercept=false)
        σᶜ = (Xₑ - Z₂*b)' * (Xₑ - Z₂*b)

        F = ((σᶜ - σ¹) * (size(Y, 1) - size(Z₂, 2))) / (size(Z, 2) * σ¹)
    end

    # return the point estimates, variance-covariance matrix,
    # and, optionally, the F-statistic from the first stage
    if Ftest
        return β, V, F
    else
        return β, V
    end

end # 2SLS()
