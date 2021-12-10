#==
    This file defines a logit regression function and all of its
    dependencies for edgel_ps2.jl
==#

using Optim

# Calculate log-likelihood at β
function likelihood(β, Y, X; intercept = true)

    # if indicated, add an intercept to the model
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    return sum(Y .* log.(exp.(X*β) ./ (1 .+ exp.(X*β))) +
                (1 .- Y) .* log.(1 ./ (1 .+ exp.(X*β))))
end # log-likelihood function

# Calculate the Hessian matrix given β
function Hessian(X, β)

    X = [ones(size(X, 1), 1) X] # add constant to X

    A = (exp.(X * β) ./ (1 .+ exp.(X * β))) .*
        (1 ./ (1 .+ exp.(X * β)))

    H = 0
    for i = 1:size(X, 1)
        H = H .+ (A[i] .* X[i, :] * transpose(X[i, :]))
    end

    return -H
end # Hessian matrix

# Define the Newton convergence algorithm
function NewtonAlg(Y, X; β₀::Matrix{Float64} = nothing,
    err::Float64 = 100.0, tol::Float64 = 1e-32, sk::Float64 = 1e-7,
    intercept::Bool = true)

    # if no initial guess is given, use all zeros
    if (is_nothing(β₀)); β₀ = zeros(size(X, 2), 1); end

    # if indicated, add an intercept to the model
    if (intercept); X = [ones(size(X, 1), 1) X]; end

    β_out = 0
    iter = 1
    while err > tol

        # update β
        β_out = β₀ - sk * inv(Hessian(X, β₀)) * transpose(score(β₀, Y, X))

        #If you have made β_out NaN, you've gone too far
        while isnan((transpose(β_out)ones(size(X, 2) + 1, 1))[1])
            sk = sk / 10
            β_out = β₀ - sk * inv(Hessian(X, β₀)) * transpose(score(β₀, Y, X))
        end

        # calculate error and update β₀
        err_new = maximum(abs.(β_out - β₀))
        β₀ = copy(β_out)
        if iter % 5 == 0
            println("Newton Iteration $(iter) with error $(err_new)")
        end
        iter += 1

        #Update sk depending on whether things are going well or not
        if err_new < err
            sk = sk * 2
        else
            sk = sk / 10
        end
        err = copy(err_new)
    end # err > tol loop

    # return converged β
    return β_out
end # Newton's algorithm

# wrapper function: logit regression
function logit(Y, X; β₀ = nothing, intercept = true, method = "BFGS", tol = 1e-32)

    # if no initial guess is given, use all zeros
    if (β₀ === nothing);
        if (intercept)
            β₀ = zeros(size(X, 2) + 1, 1)
        else
            β₀ = zeros(size(X, 2), 1)
        end
    end

    # if Newton method, use Newton's algorithm; otherwise, pass
    # method to Optim.minimize
    if method == "Newton"
        return NewtonAlg(Y, X; intercept = intercept, tol = tol)
    else
        return optimize(b -> -likelihood(b, Y, X; intercept = intercept),
            β₀, method = BFGS(), f_tol = tol, g_tol = tol).minimizer
    end

end # logit function