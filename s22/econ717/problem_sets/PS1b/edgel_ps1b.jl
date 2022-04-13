#==
    This file performs all necessary tasks for problem set 1
    of the second quarter of Econ 717

    Date created:  08 Apr 2022
    Last modified: 13 Apr 2022
    Author: Danny Edgel (edgel@wisc.edu)
==#


### Housekeeping

# load necessary packages
using Optim, Random, Plots, Printf

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

replFn = (x, i, θ) -> [j ≠ i ? θ[j] : x for j in 1:length(θ)]

for i ∈ 3:(length(θ)-1)
    ygrid = [OccChoiceObj(replFn(x, i, θ₀)) for x in 0.1:0.1:50]
    θ₀[i] = (0.1:0.1:50)[findmin(ygrid)[2]]
end

ygrid = [OccChoiceObj(replFn(x, length(θ₀), θ₀)) for x in 0.001:0.001:0.999]
θ₀[length(θ₀)] = (0.001:0.001:0.999)[findmin(ygrid)[2]]

OccChoiceObj(θ₀)
mean(SimulateData(θ₀, 1000).j .== 1)



## e) Estimate μ₁ and ρ from the simulated data, holding other params
##   fixed at their "true" values

# generate a data set with the true parameters
dat = SimulateData(θ₀, 1000)

SMMObjFun(θ, dat)

# define secondary objective function that holds other parameters fixed
obj = x -> SMMObjFun([θ₀[1], θ₀[2], x[1], θ₀[4], θ₀[5], θ₀[6], x[2]], dat; N=1000, k=2)

# derive the parameter estimate
x̂ = Optim.optimize(obj, [0.01, 0.0], [Inf, 0.99], [0.5, 0.5], Fminbox()).minimizer

θ = [θ₀[1], θ₀[2], x̂[1], θ₀[4], θ₀[5], θ₀[6], x̂[2]]

## f) Create some figures that show that x̂ identifies the true parameters

μgrid = 0.001:0.001:0.99
ρgrid = 0:0.0005:0.9995
ygridμ = [obj([x, x̂[2]]) for x in μgrid]
ygridρ = [obj([x̂[1], x]) for x in ρgrid]

plot(μgrid, ygridμ)
plot(ρgrid, ygridρ)
plot(ρgrid[ρgrid.<=0.9], ygridρ[ρgrid.<=0.9])

## g) report average number of agents who choose 1, and the
##    average and standard deviation of wages in each occupation

# simulate data with the estimated parameters
sim = SimulateData(θ, 1000)

str = @sprintf "\\begin{tabular}{r|cc} 
                    & Data & Model \\\\\\hline
    %% Choosing j=1 & %0.2f & %0.2f \\\\
    \$\\E{W|j=1}\$   & %0.2f & %0.2f \\\\
    \$\\E{W|j=2}\$   & %0.2f & %0.2f \\\\
    \$Var(W|j=1)\$   & %0.2f & %0.2f \\\\
    \$Var(W|j=2)\$   & %0.2f & %0.2f \\\\
    \end{tabular}
" mean(dat.j .== 1) mean(sim.j .== 1) mean(dat.W[dat.j.==1]) mean(sim.W[sim.j.==1]) mean(dat.W[dat.j.==2]) mean(sim.W[sim.j.==2]) var(dat.W[dat.j.==1]) var(sim.W[sim.j.==1]) var(dat.W[dat.j.==2]) var(sim.W[sim.j.==2])


open("table1.tex", "w") do io
    write(io, str)
end

## h) find w̲ such that 70% of agents choose occupation 1

obj = x -> OccChoiceObj(θ; p=0.7, w̲ = x)
ŵ̲ = (0.0:0.001:2.0)[findmin([obj(x) for x in 0.0:0.001:2.0])[2]]