#==
    This file completes all of the tasks for problem
    set 2 of the second quarter of Econ 715

    Date created:  07 Dec 2021
    Last modified: 10 Dec 2021
    Author: Danny Edgel
==#


using Printf, DataFrames, StatFiles, Optim, Distributions, Parameters, Plots

# load functions
include("functions.jl");
include("logit.jl");

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

βₓ, Vₓ = CATE(Y[subset], X[subset, 2:end], T; SE = true);

# 2a: plot CATE against experience, including 95% c.i.
Xlevels = sort(unique(X[subset, 2], dims = 1));
plot(Xlevels, βₓ, lab = "Point Estimate")
plot!(Xlevels, βₓ .+ 1.96 * Vₓ,
        fillrange = βₓ .- 1.96 * Vₓ,
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
        τ += (sum(X[subset, 2] .== Xlevels[i]) / sum(subset)) * βₓ[i]
end

# write ATE to .tex file
fname = "./econ715/Q2/problem_sets/PS2/2b.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" τ
        write(io, str)
end;


#==
        Question 3
==#

# calculate naive ATE for a single sample of n = 400
τ = simATE(Y[subset], T; B = 1, n = 400);

# write sample ATE to .tex file
fname = "./econ715/Q2/problem_sets/PS2/3.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" τ[1]
        write(io, str)
end;

# repeat for B = 500
τ = simATE(Y[subset], T; B = 500, n = 400);

# plot a histogram of ATE estimates
histogram(τ, legend = nothing, normalize = :probability)
savefig("econ715/Q2/problem_sets/PS2/fig2.png")


#==
        Question 4
==#

# calculate ATE for a single sample of n = 400
τₓ, τ = simCATE(Y[subset], X[subset, 2:end], T; B = 1, n = 400);

# write sample ATE to .tex file
fname = "./econ715/Q2/problem_sets/PS2/4.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" τ[1]
        write(io, str)
end;

# repeat for B = 500
τₓ, τ = simCATE(Y[subset], X[subset, 2:end], T; B = 500, n = 400);

# plot a histogram of ATE estimates
histogram(τ, legend = nothing, normalize = :probability)
savefig("econ715/Q2/problem_sets/PS2/fig3.png")


#==
        Question 5
==#

# estimate a propensity score for the treatment using logit 
τ = simPS(Y[subset], X[subset, 2:3], T; n = 400, B = 1);

# write sample ATE to .tex file
fname = "./econ715/Q2/problem_sets/PS2/5.tex";
open(fname, "w") do io
        str = @sprintf "%1.3f" τ[1]
        write(io, str)
end;

# repeat for B = 500
τ = simPS(Y[subset], X[subset, 2:3], T; n = 400, B = 500);

# plot a histogram of ATE estimates
histogram(τ, legend = nothing, normalize = :probability)
savefig("econ715/Q2/problem_sets/PS2/fig4.png")