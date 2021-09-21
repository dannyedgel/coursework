using Parameters, Plots, LinearAlgebra, StatsPlots #import the libraries we want
cd("C:/Users/edgel/Google Drive/UW-Madison/f21/econ899/Q1/problem_sets/PS2")
include("edgel_model_functions.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func, μ = res
@unpack a_grid = prim

##############Make plots

#stationary distribution 
Q = sum(res.μ, dims = 1)
Plots.plot(a_grid, transpose(Q))

#value function
Plots.plot(k_grid, transpose(val_func), title="Value Functions",
                label = ["S = e" "S = u"], legend=:bottomleft)
Plots.savefig("02_Value_Functions.png")

#policy functions
Plots.plot(k_grid, transpose(pol_func), title="Policy Functions",
                label = ["Z = 1.25" "Z = 0.2"], legend=:topleft)
Plots.savefig("02_Policy_Functions.png")

#changes in policy function
pol_func_δ = transpose(copy(pol_func)[:,2:end].-copy(pol_func)[:,1:(end-1)])
Plots.plot(k_grid[1:(end-1)], pol_func_δ, title="Policy Functions Changes",
                label = ["Z = 1.25" "Z = 0.2"], legend=:bottomright)
Plots.savefig("02_Policy_Functions_Changes.png")

println("All done!")
################################
