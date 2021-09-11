using Parameters, Plots #import the libraries we want
include("edgel_growth_model.jl") #import the functions that solve our growth model
cd("C:/Users/edgel/Google Drive/UW-Madison/f21/econ899/Q1/problem_sets/PS1")

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func = res
@unpack k_grid = prim

##############Make plots
#value function
Plots.plot(k_grid, transpose(val_func), title="Value Functions",
                label = ["Z = 1.25" "Z = 0.2"], legend=:bottomleft)
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
