#==
    This file performs all necessary tasks for problem set 1
    of the second quarter of Econ 717

    Date created:  08 Apr 2022
    Last modified: 08 Apr 2022
    Author: Danny Edgel (edgel@wisc.edu)
==#


### Housekeeping

# load necessary packages
using Optim, Random

# declare user-written functions 
include("./ps1b_functions.jl")

# set seed for reproducibility
Random.seed!(115)

### Questions

## b) Find θ such that ~60% of simulated agents choose 1

# initialize the model with default parameters
prim, θ = Initialize()

# define bounds for search
θ_lb = [0, 0, -Inf, -Inf, 0, 0, 0]
θ_ub = [Inf, Inf, Inf, Inf, Inf, Inf, 0.99]

# beginning with θ, perform a grid search over the parameter
# space to find θ such that ~60% of simulated agents choose 1
θ₀ = Optim.optimize(OccChoiceObj, θ_lb, θ_ub, θ, Fminbox()).minimizer
