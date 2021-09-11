using Parameters, Plots #import the libraries we want
include("edgel_growth_model.jl") #import the functions that solve our growth model
cd("C:/Users/edgel/Google Drive/UW-Madison/f21/econ899/Q1/problem_sets/PS1")

prim, res = Initialize() #initialize primitive and results structs
@elapsed Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func = res
@unpack k_grid = prim

### NEXT: CHANGE PLOTS

##############Make plots
#value function
Plots.plot(k_grid, val_func, title="Value Function")
Plots.savefig("02_Value_Functions.png")

#policy functions
Plots.plot(k_grid, pol_func, title="Policy Functions")
Plots.savefig("02_Policy_Functions.png")

#changes in policy function
pol_func_δ = copy(pol_func).-k_grid
Plots.plot(k_grid, pol_func_δ, title="Policy Functions Changes")
Plots.savefig("02_Policy_Functions_Changes.png")

println("All done!")
################################
