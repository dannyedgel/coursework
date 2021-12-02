#==
    This file defines all of the functions used by edgel_ps1.jl
==#
using LinearAlgebra, Statistics

# function that simply runs OLS and outputs coefficients and
# variance-covariance matrix
function OLS(Y, X; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # calculate coefficients 
    β = inv(X' * X) * X' * Y

    # calculate variance-covariance matrix
    Ω = (Y - X * β)' * (Y - X * β) / (size(Y, 1) - size(X, 2))
    V = (X' * X) .* Ω

    # return coefficients and variance-covariance matrix
    return β, V

end # OLS function 

# criterion function for quantile regression
function QuantCriterion(Y, X, β, τ; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # Calculate ρ function for each observation
    ρ = (τ - 1) * (Y - X * β) .* (Y .≤ X * β) + τ * (Y - X * β) .* (Y .> X * β)

    # Output mean value
    return (1 / size(Y, 1)) * sum(ρ)

end # Quantile criterion function

# function that generates bootstrap estimates for a quantile 
# regression coefficient
function QBootstap(Y, X, β; B::Int64 = 1000, τ = 0.75, intercept = true)

    # generate bootstrap samples of Y and X
    Ys, Xs = Y[rand(1:length(Y), (length(Y), B))], X[rand(1:length(Y), (length(Y), B)), :]

    # initialize vector of θhats
    Θhat = zeros(length(β), B)

    # generate vector of θhats
    for b = 1:B
        bhat = optimize(k -> QuantCriterion(Ys[:, b], Xs[:, b, :], k, τ, intercept = intercept), β, LBFGS()).minimizer
        Θhat[:, b] = bhat
    end

    # return bootstrap SE
    return ((B - 1)^(-1)) * sum((Θhat .- mean(Θhat, dims = 2)) .^ 2, dims = 2)

end # Bootstrap SE

# GMM objective function for quantile regression
function Qobj(Y, X, β, τ)

    # define the objective function for a given β
    W = Matrix{Float64}(I, size(X, 2), size(X, 2))

    # initialize the sum for Γ
    Γ = zeros(size(X[1, :]))
    for i = 1:size(X, 2)
        Γ .+= X[i, :] * ((Y[i] ≤ dot(X[i, :], β)) - τ)
    end
    Γ = (1 / length(Y)) * Γ

    # Output mean value
    return transpose(Γ) * W * Γ

end # Qobj function

# function that runs a quantile regression with GMM
function QGMM(Y, X; τ = .5, intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # run the optimization
    β, V = OLS(Y, X; intercept = false)
    βhat = optimize(b -> Qobj(Y, X, b, τ), β, LBFGS()).minimizer

    # estimate the variance-covariance matrix
    ε = Y - X * βhat
    disthat = fit(Normal, ε)
    Γ, Ω = zeros(size(X, 2), size(X, 2)), zeros(size(X, 2), size(X, 2))
    for i = 1:length(Y)
        Γ .+= X[i, :] * transpose(X[i, :]) * pdf(disthat, 0)
        Ω .+= X[i, :] * transpose(X[i, :])
    end
    Γ = (1 / length(Y)) * Γ
    Ω = (1 / length(Y)) * Ω * τ * (1 - τ)


    # return the coefficient and variance-covariance matrix
    return βhat, inv(transpose(Γ) * inv(Ω) * Γ)

end # GMM quantile regression 


# function that estimates the model with FGLS
function FGLS(Y, X; τ = 0.5, intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # generate array of data for each unique education-experience pair
    # NOTE: ω omitted because I only included it out of ignorance
    uqX = unique(X, dims = 1)
    CellDat = Array{Any}(zeros(size(uqX, 1), 1))
    α = Array{Any}(zeros(size(uqX, 1), 1))
    ω = Array{Any}(zeros(size(uqX, 1), 1))
    Yhat = Array{Float64}(zeros(size(uqX, 1), 1))
    exclude = Array{Int64}(undef, 0)
    for i = 1:size(uqX, 1)
        CellDat[i] = Y[sum(X .== uqX[i, :]', dims = 2).==size(X, 2), :]
        if length(CellDat[i]) > 2
            α[i], ω[i] = QGMM(CellDat[i], ones(size(CellDat[i])); τ = τ, intercept = false) 
            ω[i] = ω[i][1]
            Yhat[i] = mean(CellDat[i]) + ω[i] * α[i][1]
        else
            push!(exclude, i)
        end
    end

    # subset to only include those that were not excluded
    if length(exclude) > 0
        Yhat = Yhat[1:length(Yhat).∉Ref(Set(exclude))]
        uqX = uqX[1:size(uqX, 1).∉Ref(Set(exclude)), :]
        ω = ω[1:length(ω).∉Ref(Set(exclude))]
    end

    # generate Ω matrix for the full sample
    ω = convert(Array{Float64}, ω)
    Ω = ω.*Matrix(I, length(ω), length(ω))

    # calculate FGLS coefficients and variance-covariance matrix
    β = inv(uqX' * inv(Ω) * uqX) * uqX' * inv(Ω) * Yhat
    ε = Yhat - uqX * β
    V = inv(uqX' * inv(Ω) * uqX) * (uqX' * inv(Ω) * ε * ε' * inv(Ω) * uqX) * inv(uqX' * inv(Ω) * uqX)

    # return coefficients and variance-covariance matrix
    return β, V

end # FGLS function

# Object that stores CI arrays from sample simulations
struct CI
    QBoot::Array{Float64,2}
    GMM::Array{Float64,2}
    FGLS::Array{Float64,2}
end

# function that simulates a subsample of the data and
# computes confidence intervals for the coefficients β₁^.75
function SimCI(Y, X; τ = 0.75, n = 400, J = 1000, intercept = true)

    # if indicated, add an intercept to the model
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    # generate subsamples of the data
    Ys = Y[rand(1:length(Y), (n, J))]
    Xs = X[rand(1:length(Y), (n, J)), :]

    # initialize array of confidence intervals for each estimation method
    SampleCIs = CI(zeros(J, 2), zeros(J, 2), zeros(J, 2))

    # calculate CIs for each estimation method
    β, V = OLS(Y, X, intercept = false) # starting value for β
    for j = 1:J
    
        # Quantile regression with bootstrap
        βQ = optimize(b -> QuantCriterion(Ys[:, j], Xs[:, j, :], b, τ; intercept = false), β, LBFGS()).minimizer
        QBstr = QBootstap(Ys[:, j], Xs[:, j, :], β; B = 10, τ = τ, intercept = false)
    
        SampleCIs.QBoot[j, 1] = βQ[2] - QBstr[2]
        SampleCIs.QBoot[j, 2] = QBstr[2] - βQ[2]
    
        # GMM 
        βG, VG = QGMM(Ys[:, j], Xs[:, j, :]; τ = τ, intercept = false)
        SampleCIs.GMM[j, 1] = βG[2] - sqrt(VG[2, 2])
        SampleCIs.GMM[j, 2] = sqrt(VG[2, 2]) + βG[2]
    
        # FGLS
        βF, VF = FGLS(Ys[:, j], Xs[:, j, :]; τ = τ, intercept = false)
        SampleCIs.FGLS[j, 1] = βF[2] - sqrt(VF[2, 2])
        SampleCIs.FGLS[j, 2] = sqrt(VF[2, 2]) + βF[2]
    
    end


    # return CI struct 
    return SampleCIs
end # SimCI function