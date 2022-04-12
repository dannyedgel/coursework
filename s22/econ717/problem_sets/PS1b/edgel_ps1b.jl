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
θ_lb = [0, 0, 0, 0, 0, 0, 0]
θ_ub = [Inf, Inf, Inf, Inf, Inf, Inf, 0.99]


# beginning with θ, perform a grid search over the parameter
# space (parameter-by-parameter) to find θ such that ~60% of 
# simulated agents choose 1

θ₀ = θ

replFn = (x, i, θ) -> [j≠i ? θ[j] : x for j in 1:length(θ)]

for i∈3:(length(θ) - 1)
    ygrid = [OccChoiceObj(replFn(x, i, θ₀)) for x in 0.1:0.1:50]
    θ₀[i] = (0.1:0.1:50)[findmin(ygrid)[2]]
end

ygrid = [OccChoiceObj(replFn(x, length(θ₀), θ₀)) for x in 0.001:0.001:.999]
θ₀[length(θ₀)] = (0.001:0.001:.999)[findmin(ygrid)[2]]

OccChoiceObj(θ₀)
mean(SimulateData(θ₀, 1000).j .== 1)



## e) Estimate μ₁ and ρ from the simulated data, holding other params
##   fixed at their "true" values

# generate a data set with the true parameters
dat = SimulateData(θ₀, 1000)

SMMObjFun(θ, dat)

# define secondary objective function that holds other parameters fixed
obj = x -> SMMObjFun([θ₀[1], θ₀[2], x[1], θ₀[4], θ₀[5], θ₀[6], x[2]], dat)

# derive the parameter estimate
x̂ = Optim.optimize(obj, [0.01, 0.0], [Inf, .99], [.5, .5], Fminbox()).minimizer


## f) Create some figures that show that x̂ identifies the true parameters


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

