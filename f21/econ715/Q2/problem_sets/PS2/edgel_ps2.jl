#==
    This file completes all of the tasks for problem
    set 2 of the second quarter of Econ 715

    Date created:  07 Dec 2021
    Last modified: 09 Dec 2021
    Author: Danny Edgel
==#


using Printf, DataFrames, StatFiles, Optim, Distributions, Parameters, Plots

# load functions
include("functions.jl")

# Load 2009 CPS from Bruce's website
df = DataFrame(StatFiles.load("econ715/Q2/problem_sets/PS2/data/cps09mar.dta"))

# subset and mutate as directed
df.exp = df.age - df.education .- 6;
df = df[(df.exp.>=9).&(df.exp.<=14).&(df.education.>=11).&(df.education.<=17), :]


#==
        Question 1
==#

# regress log(earnings) on education, experience, and exp^2
Y = log.(df[!, [:earnings]]) |> Matrix
X = df[!, [:education, :exp,]] |> Matrix;
X = identity.([X (X[:, 2] .^ 2)]);

β, V = OLS(Y, X);

# write education coefficient to .tex file
fname = "./econ715/Q2/problem_sets/PS2/1.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" β[2]
        write(io, str)
end;


#==
        Question 2
==#

# subset the data to observations with education = 12 or = 16
subset = (X[:, 1] .== 12) .| (X[:, 1] .== 16);
#Y = Y[subset]; X = X[subset, :];

# generate treatment treatment variable
T = X[subset, 1] .== 16;

# run OLS on the treatment interacted with experience to obstain the CATE
B = [T (X[subset, 2] .* T) (X[subset, 3] .* T) X[subset, 2] X[subset, 3]]
βₜ, Vₜ = OLS(Y[subset], B);

# for each value of experience, calculate CATE
Xlevels = sort(unique(X[subset, 2]))
CATE = βₜ[2] .+ βₜ[3] * Xlevels + βₜ[4] *(Xlevels^2);

# calculate SE of each CATE 
seCATE = sqrt.(Vₜ[2, 2] .+ Vₜ[3, 3] * (Xlevels .^ 2) .+ 2 * Xlevels * Vₜ[2, 3])

# 2a: plot CATE against experience, including 95% c.i.
plot(Xlevels, CATE, lab = "Point Estimate")
plot!(Xlevels, CATE .+ 1.96 * seCATE,
        fillrange = CATE .- 1.96 * seCATE,
        fillalpha = 0.35,
        lw = 0,
        fillcolor = :grey,
        linecolor = nothing,
        xlabel = "Experience (years)",
        ylabel = "log(Wage)",
        #title = "CATE of increasing education from 12 years to 16",
        lab = "95% c.i.")
savefig("econ715/Q2/problem_sets/PS2/fig1.png")

# 2b: estimate ATE
τ = 0;
for i ∈ 1:length(Xlevels)
        τ += (sum(X[subset, 2] .== Xlevels[i]) / sum(subset)) * CATE[i]
end

# write ATE to .tex file
fname = "./econ715/Q2/problem_sets/PS2/2b.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" τ
        write(io, str)
end;