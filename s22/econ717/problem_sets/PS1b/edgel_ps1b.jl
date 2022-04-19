#==
    This file performs all necessary tasks for problem set 1
    of the second quarter of Econ 717

    Date created:  08 Apr 2022
    Last modified: 19 Apr 2022
    Author: Danny Edgel (edgel@wisc.edu)
==#


### Housekeeping

# load necessary packages
using Optim, Random, Plots, Printf

# declare user-written functions 
include("./ps1b_functions.jl")

# choose number of simulations to run
n = 10000

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

θ₀ = copy(θ)

replFn = (x, i, θ) -> [j ≠ i ? θ[j] : x for j in 1:length(θ)]
range = 0.1:0.005:10
iter = 0
while (abs(mean(SimulateData(θ₀, n; seed=3609).j .== 1) - 0.6) > 0.005) & (iter < 5)
    for i ∈ 3:(length(θ)-1)
        ygrid = [OccChoiceObj(replFn(x, i, θ₀)) for x in range]
        θ₀[i] = (range)[findmin(ygrid)[2]]
    end

    ygrid = [OccChoiceObj(replFn(x, length(θ₀), θ₀)) for x in -0.999:0.001:0.999]
    θ₀[length(θ₀)] = (-0.999:0.001:0.999)[findmin(ygrid)[2]]

    iter += 1
    print("\nIteration " * string(iter) * ", err = " *
          string(abs(mean(SimulateData(θ₀, n; seed=3609).j .== 1) - 0.6)))
end

OccChoiceObj(θ₀)
mean(SimulateData(θ₀, n; seed=3609).j .== 1)

str = @sprintf "\$\\theta=(1, 1, %.2f, %.2f, %.2f, %.2f, %.2f)\$ results in \$%.2f\$\\%% of simulated agents choosing 1 for the chosen randomization seed.
" θ₀[3] θ₀[4] θ₀[5] θ₀[6] θ₀[7] mean(SimulateData(θ₀, n; seed=3609).j .== 1) * 100

open("c.tex", "w") do io
    write(io, str)
end

## e) Estimate μ₁ and ρ from the simulated data, holding other params
##   fixed at their "true" values

# generate a data set with the true parameters
dat = SimulateData(θ₀, n)

SMMObjFun(θ, dat)

# define secondary objective function that holds other parameters fixed
obj = x -> SMMObjFun([θ₀[1], θ₀[2], x[1], θ₀[4], θ₀[5], θ₀[6], x[2]], dat;
    N=Int(n), m=5)

# derive the parameter estimate
x̂ = Optim.optimize(obj, [0.01, -0.99], [Inf, 0.99], [θ[3], θ[7]], Fminbox(LBFGS())).minimizer

str = @sprintf "\$(\\hat{\\mu_1}, \\hat{\\rho})=(%.2f, %.2f)\$.
" x̂[1] x̂[2]

open("e.tex", "w") do io
    write(io, str)
end

## f) Create some figures that show that x̂ identifies the true parameters

μgrid = 0.001:0.001:0.99
ρgrid = -0.1:0.0005:0.5
#ygridμ = [obj([x, x̂[2]]) for x in μgrid]
#ygridρ = [obj([x̂[1], x]) for x in ρgrid]
ygridμ = [obj([x, θ₀[7]]) for x in μgrid]
ygridρ = [obj([θ₀[3], x]) for x in ρgrid]

plot(μgrid, ygridμ, labels="SMM Obj Fun", title="\\mu_1 identification")
plot!([θ₀[3]], seriestype="vline", labels="True")
plot!([x̂[1]], seriestype="vline", labels="Estimate", linestyle=:dash)
savefig("f1.png")

plot(ρgrid, ygridρ, labels="SMM Obj Fun", title="\\rho identification")
plot!([θ₀[7]], seriestype="vline", labels="True")
plot!([x̂[2]], seriestype="vline", labels="Estimate", linestyle=:dash)
savefig("f2.png")


θ = [θ₀[1], θ₀[2], x̂[1], θ₀[4], θ₀[5], θ₀[6], x̂[2]]


## g) report average number of agents who choose 1, and the
##    average and standard deviation of wages in each occupation

# simulate data with the estimated parameters (table defined in question h)
sim = SimulateData(θ, n; seed=1326)


## h) find w̲ such that 70% of agents choose occupation 1
obj = x -> OccChoiceObj(θ; p=0.7, w̲=x, N=Int(n / 100), seed=214)
ŵ̲ = (0.0:0.001:2.0)[findmin([obj(x) for x in 0.0:0.001:2.0])[2]]

str = @sprintf "\$%.3f\$." ŵ̲
open("h.tex", "w") do io; write(io, str); end


sim2 = SimulateData(θ, n; seed=1326, w̲=ŵ̲)

str = @sprintf "\\begin{tabular}{r|ccc} 
                        & Data  & Model & \$\\underline{\\hat{W}}_1\$ \\\\\\hline
    \\%% Choosing j=1   & %0.3f & %0.3f & %0.3f                       \\\\
    \$\\E{W|j=1}\$      & %0.3f & %0.3f & %0.3f                       \\\\
    \$\\E{W|j=2}\$      & %0.3f & %0.3f & %0.3f                       \\\\
    \$Var(W|j=1)\$      & %0.3f & %0.3f & %0.3f                       \\\\
    \$Var(W|j=2)\$      & %0.3f & %0.3f & %0.3f                       \\\\
    \\end{tabular}
" mean(dat.j .== 1) mean(sim.j .== 1) mean(sim2.j .== 1) mean(dat.W[dat.j.==1]) mean(sim.W[sim.j.==1]) mean(sim2.W[sim2.j.==1]) mean(dat.W[dat.j.==2]) mean(sim.W[sim.j.==2]) mean(sim2.W[sim2.j.==2]) var(dat.W[dat.j.==1]) var(sim.W[sim.j.==1]) var(sim2.W[sim2.j.==1]) var(dat.W[dat.j.==2]) var(sim.W[sim.j.==2]) var(sim2.W[sim2.j.==2])

open("table1.tex", "w") do io
    write(io, str)
end