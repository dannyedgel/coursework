#==
    This file defines all of the functions used by edgel_ps1.jl
==#
using LinearAlgebra, Statistics

# function that simply runs OLS and outputs coefficients and
# variance-covariance matrix
function OLS(Y, X; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    # calculate coefficients 
    β = inv(X' * X) * X' * Y;

    # calculate variance-covariance matrix
    Ω = (Y - X*β)'*(Y - X*β) / (size(Y, 1) - size(X, 2));
    V = (X'*X) .* Ω;

    # return coefficients and variance-covariance matrix
    return β, V

end # OLS function 

# criterion function for quantile regression
function QuantCriterion(Y, X, β, τ; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    # Calculate ρ function for each observation
    ρ = (τ - 1)*(Y - X*β).*(Y .≤ X*β) + τ*(Y - X*β).*(Y .> X*β)

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
        bhat = optimize(k -> QuantCriterion(Ys[:, b], Xs[:, b, :], k, τ), β, LBFGS()).minimizer
        Θhat[:, b] = bhat
    end

    # return bootstrap SE
    return ((B - 1)^(-1)) * sum((Θhat .- mean(Θhat, dims = 2)).^2, dims = 2)

end # Bootstrap SE

# GMM objective function for quantile regression
function Qobj(Y, X, β, τ)

    # define the objective function for a given β
    W = Matrix{Float64}(I, size(X, 2), size(X, 2)) 

    # initialize the sum for Γ
    Γ = zeros(size(X[1, :]));
    for i = 1:size(X, 2)
        Γ .+= X[i, :]*((Y[i] ≤ dot(X[i, :], β)) - τ)
    end
    Γ = (1/length(Y))*Γ

    # Output mean value
    return transpose(Γ) * W * Γ

end # Qobj function

# function that runs a quantile regression with GMM
function QGMM(Y, X, τ; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # run the optimization
    β, V = OLS(Y, X, intercept = false)
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