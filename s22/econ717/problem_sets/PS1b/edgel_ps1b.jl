#==
    This file performs all necessary tasks for problem set 1
    of the second quarter of Econ 717

    Date created:  08 Apr 2022
    Last modified: 10 Apr 2022
    Author: Danny Edgel (edgel@wisc.edu)
==#


### Housekeeping

# load necessary packages
using Optim, Random

# declare user-written functions 
include("./ps1b_functions.jl")

### Questions

## b) See SimulateData() in ps1b_functions.jl

## c) Find θ such that ~60% of simulated agents choose 1

# initialize the model with default parameters
prim, θ = Initialize()

# define bounds for search
θ_lb = [0, 0, -Inf, -Inf, 0, 0, 0]
θ_ub = [Inf, Inf, Inf, Inf, Inf, Inf, 0.99]


# beginning with θ, perform a grid search over the parameter
# space to find θ such that ~60% of simulated agents choose 1

θ₀ = Optim.optimize(OccChoiceObj, θ_lb, θ_ub, θ, Fminbox()).minimizer

## e) Estimate μ₁ and ρ from the simulated data, holding other params
##   fixed at their "true" values




#==


test = x -> OccChoiceObj([1, 1, x[1], x[2], x[3], x[4], x[5]])


θ₀ = Optim.optimize(test, θ_lb[3:end], θ_ub[3:end], θ[3:end], Fminbox())


bounds = [(θ_lb[i], θ_ub[i]) for i∈1:length(θ)]

res = BlackBoxOptim.bboptimize(OccChoiceObj, Method = :random_search,
                                SearchRange=bounds, NumDimensions=length(θ))
best_candidate(res)

opt = Opt(:LD_MMA, length(θ))
opt.lower_bounds, opt.upper_bounds  = θ_lb, θ_ub
opt.min_objective = OccChoiceObj

(minf,minx,ret) = NLopt.optimize!(opt, θ)
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")

