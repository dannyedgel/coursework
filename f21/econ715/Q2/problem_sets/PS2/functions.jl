#==
    This file defines all of the functions used by edgel_ps2.jl
==#
using LinearAlgebra, Statistics, StatsBase

# function that simply runs OLS and outputs coefficients and
# variance-covariance matrix
function OLS(Y, X; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept)
        X = [ones(size(X, 1), 1) X]
    end

    # calculate coefficients 
    β = inv(X'*X)*X'*Y

    # calculate variance-covariance matrix
    Ω = (Y - X*β)'*(Y - X*β)/(size(Y, 1) - size(X, 2))
    V = inv(X'*X).*Ω

    # return coefficients and variance-covariance matrix
    return β, V

end # OLS function 

# function that calculates the CATE for each discrete X
# and optionally returns the SE of the CATE
function CATE(Y, X, T; SE = false)

    # generate new matrix with interactions between T and each
    # column of X
    Xₜ = T
    if size(X, 2) > 1
        for i = 1:size(X, 2)
            Xₜ = [Xₜ (X[:, i] .* T)]
        end
    end
    for i = 1:size(X, 2)
        Xₜ = [Xₜ X[:, i]]
    end

    # run OLS on the new matrix
    βₜ, Vₜ = OLS(Y, Xₜ)

    # calculate the CATE
    Xlevels = unique(X, dims = 1); 
    Xlevels = Xlevels[sortperm(Xlevels[:, 1]), :]
    τₓ = βₜ[2] .+ Xlevels*βₜ[3:(size(X, 2)+2)]

    # if indicated, calculate the standard error
    if (SE)
        W = [zeros(size(Xlevels, 1), 1) ones(size(Xlevels, 1), 1) Xlevels]
        W = [W zeros(size(Xlevels, 1), length(βₜ) - size(W, 2))];
        Ωₓ = W*Vₜ*W'; SEₓ = sqrt.(diag(Ωₓ));

        # return the CATE and standard error
        return τₓ, SEₓ
    else
        # return the CATE
        return τₓ
    end

end # CATE function


# function that calculates the naive ATE for subsamples of the data, B times 
function simATE(Y, T; B = 500, n = 400, SE = false, seed = 4231)

    # sample n observations without replacement, B times
    Random.seed!(seed)
    S = Array{Int64}(zeros(n, B));
    for i = 1:B; S[:, i] = sample(1:length(Y), n, replace = false); end

    # initalize array to hold ATEs (and, optionally, SEs)
    if (SE)
        τ = Array{Float64}(zeros(B, 1));
        V = Array{Float64}(zeros(B, 1));
    else
        τ = Array{Float64}(zeros(B, 1));
    end

    # loop through simulations, calculating ATE and SE
    for i = 1:B
        # run OLS on the subsample
        βᵢ, Vᵢ = OLS(Y[S[:, i]], T[S[:, i]])

        # calculate ATE and SE
        if (SE)
            τ[i], V[i] = βᵢ[2], sqrt(Vᵢ[2, 2])
        else
            τ[i] = βᵢ[2]
        end
    end

    # return results
    if (SE); return τ, V; else; return τ; end

end # simATE function

# function that calculates the CATE for subsamples of the data, B times
function simCATE(Y, X, T; B = 500, n = 400, SE = false, seed = 4231)

    # define simple propensity score function
    pscore =

    # sample n observations without replacement, B times
    Random.seed!(seed)
    S = Array{Int64}(zeros(n, B))
    for i = 1:B; S[:, i] = sample(1:length(Y), n, replace = false); end

    # initialize array to hold CATEs and ATEs (and, optionally, SEs)
    Xlevels = unique(X, dims = 1)
    if (SE)
        Θ = Array{Float64}(zeros(size(Xlevels, 1), B))
        τ = Array{Float64}(zeros(B, 1))
        Ω = Array{Float64}(zeros(size(Xlevels, 1), B))
    else
        Θ = Array{Float64}(zeros(size(Xlevels, 1), B))
        τ = Array{Float64}(zeros(B, 1))
    end

    # loop through simulations, calculating CATEs
    for i = 1:B
        # calculate CATEs for each subsample
        if (SE)
            Θ[:, i], Ω[:, i] = CATE(Y[S[:, i]], X[S[:, i], :],
                T[S[:, i]], SE = true)
        else
            Θ[:, i] = CATE(Y[S[:, i]], X[S[:, i], :],
                T[S[:, i]], SE = false)
        end

        # calculate ATEs for each subsample
        τ[i] = 0
        for j = 1:size(Xlevels, 1)
            τ[i] += (sum(X[S[:, i], :] .== Xlevels[j]) / n) * Θ[j, i]
        end
    end

    # return CATEs and ATEs (and, optionally, SEs)
    if (SE)
        return Θ, τ, Ω
    else
        return Θ, τ
    end

end # simCATE function

# simple propensity score calculation function
function pscore(T, X; intercept = true)

    # if indicated, add an intercept to the X array
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    # return propensity score
    return exp.(X * logit(T, X, intercept = false)) ./ (1 .+ exp.(X * logit(T, X, intercept = false)))

end # pscore function

# function that simulates the propensity score estimation of ATE
function simPS(Y, X, T; B = 500, n = 400, seed = 4231)

    # sample n observations without replacement, B times
    Random.seed!(seed)
    S = Array{Int64}(zeros(n, B));
    for i = 1:B; S[:, i] = sample(1:length(Y), n, replace = false); end

    # initialize array to hold PSs
    τ = Array{Float64}(zeros(B, 1));

    # loop through simulations, calculating PS and ATE
    for i = 1:B
        # calculate propensity score
        pₓ = pscore(T[S[:, i]], X[S[:, i], :])
    
        # sort observations on the propensity score
        Xᵢ = unique([X[S[:, i], :] pₓ], dims = 1)
        Xᵢ = Xᵢ[sortperm(Xᵢ[:, end]), :]
        nᵢ = size(Xᵢ, 1)
    
        # calculate ATE for two highest, lowest, and median propensity scores
        lowest = (X[S[:, i], :] .== Xᵢ[1, 1:(end-1)]') .| 
                    (X[S[:, i], :] .== Xᵢ[2, 1:(end-1)]')
        highest = (X[S[:, i], :] .== Xᵢ[end-1, 1:(end-1)]') .|
                    (X[S[:, i], :] .== Xᵢ[end, 1:(end-1)]')
        middle = (X[S[:, i], :] .== Xᵢ[Int(floor(nᵢ/2)), 1:(end-1)]') .|
                    (X[S[:, i], :] .== Xᵢ[Int(floor(nᵢ/2))+1, 1:(end-1)]')

        # adjust indexes to a single index vector
        lowest = (sum(lowest, dims = 2) .== size(X, 2))[:, 1]
        middle = (sum(middle, dims = 2) .== size(X, 2))[:, 1]
        highest = (sum(highest, dims = 2) .== size(X, 2))[:, 1]
    
        τₗ = simATE(Y[S[:, i]][lowest], T[S[:, i]][lowest], 
                B = 1, n = sum(lowest))[1]
        τₘ = simATE(Y[S[:, i]][middle], T[S[:, i]][middle], 
                B = 1, n = sum(middle))[1]
        τₕ = simATE(Y[S[:, i]][highest], T[S[:, i]][highest], 
                B = 1, n = sum(highest))[1]
    
        # calculate overall ATE
        τ[i] = (sum(lowest) / n) * τₗ + (sum(middle) / n) * τₘ + (sum(highest) / n) * τₕ
    end # end simulation loop

    # return results
    return τ

end # simPS function