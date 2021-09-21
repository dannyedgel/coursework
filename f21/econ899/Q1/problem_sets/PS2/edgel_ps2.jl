using Parameters, Plots, LinearAlgebra, StatsPlots #import the libraries we want
cd("C:/Users/edgel/Google Drive/UW-Madison/f21/econ899/Q1/problem_sets/PS2")
include("edgel_model_functions.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, μ, q̄ = res
@unpack a_grid = prim

##############Make plots

#stationary distribution 
Q = sum(res.μ, dims = 1)
Plots.plot(a_grid, transpose(Q))
Plots.plot(a_grid, transpose(μ))

#value function
Plots.plot(a_grid, transpose(val_func), title="Value Functions",
                label = ["S = e" "S = u"], legend=:topleft)
Plots.savefig("02_Value_Functions.png")

#policy functions
Plots.plot(a_grid, transpose(pol_func), title="Policy Functions",
                label = ["S = e" "S = u"], legend=:topleft)
Plots.savefig("02_Policy_Functions.png")

println("All done!")
################################
