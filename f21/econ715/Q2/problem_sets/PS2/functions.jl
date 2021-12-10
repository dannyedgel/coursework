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
    β = inv(X' * X) * X' * Y

    # calculate variance-covariance matrix
    Ω = (Y - X*β)'*(Y - X * β)/(size(Y, 1) - size(X, 2))
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
function simATE(Y, T; B = 500, n = 400, SE = false)
    # sample n observations without replacement, B times
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
function simCATE(Y, X, T; B = 500, n = 400, SE = false)
    # sample n observations without replacement, B times
    S = Array{Int64}(zeros(n, B));
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
    if (SE); return Θ, τ, Ω; else; return Θ, τ; end
    
end # simCATE function