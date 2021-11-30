#==
    This file completes all of the tasks for problem
    set 1 of the second quarter of Econ 715

    Date created:  18 Nov 2021
    Last modified: 29 Nov 2021
    Author: Danny Edgel
==#


using Printf, DataFrames, StatFiles, Optim, Distributions

# load functions
include("functions.jl")

# Load 2009 CPS from Bruce's website
df = DataFrame(StatFiles.load("econ715/Q2/problem_sets/PS1/data/cps09mar.dta"))

# subset and mutate as directed
df.exp = df.age - df.education .- 6;
df = df[(df.exp.>=9).&(df.exp.<=14).&(df.education.>=11).&(df.education.<=17), :]


#==
        Question 1
==#

## (a) regress log(earnings) on education, experience, and exp^2
Y = log.(df[!, [:earnings]]) |> Matrix
X = df[!, [:education, :exp,]] |> Matrix;
X = identity.([X (X[:, 2] .^ 2)]);

β, V = OLS(Y, X);

# write education coefficient to .tex file
fname = "./econ715/Q2/problem_sets/PS1/1a.tex";
open(fname, "w") do io
    str = @sprintf "%1.3f" β[2]
    write(io, str)
end;



## (b) solve for the .5 and .75 quantiles, using β from (a) as 
## the starting value in the BFGS algorithm 

# estimate coefficients
β₁ = optimize(b -> QuantCriterion(Y, X, b, 0.5), β, LBFGS()).minimizer;
β₂ = optimize(b -> QuantCriterion(Y, X, b, 0.75), β, LBFGS()).minimizer;

# write education coefficients to .tex file
fname = "./econ715/Q2/problem_sets/PS1/1b.tex";
open(fname, "w") do io
    str = @sprintf "\\bhat_1^{.5} = %1.3f, \\quad\\quad \\bhat_1^{.75} = %1.3f" β₁[2] β₂[2]
    write(io, str)
end;

## (c) calculate the bootstrap standard error for β₂
se = QBootstap(Y, X, β; B = 10);

# write bootstrap SE of education coefficient to .tex file
fname = "./econ715/Q2/problem_sets/PS1/1c.tex";
open(fname, "w") do io
    str = @sprintf "%1.8f" se[2]
    write(io, str)
end;

#==
        Question 2
==#

## run quantile regression for .5 and .75 quantiles using 
## GMM with the identity matrix
β₁, V₁ = QGMM(Y, X, 0.5);
β₂, V₂ = QGMM(Y, X, 0.75);

# write education coefficient and SE for τ=.75 to .tex file
fname = "./econ715/Q2/problem_sets/PS1/2e.tex";
open(fname, "w") do io
    str = @sprintf "\\bhat_1^{.75} = %1.3f, \\quad\\quad se(\\bhat_1^{.75}) = %1.3f" β₂[2] sqrt(V₂[2, 2])
    write(io, str)
end;


#==
        Question 3
==#

# generate array of data for each unique education-experience pair
uqX = unique(X, dims = 1)
CellDat = Array{Any}(size(uqX, 1)); α = zeros(size(uqX, 1), 1);
V = zeros(size(uqX, 1), 1);
for i = 1:size(uqX, 1)
    CellDat[i] = Y[X==uqX[i], :]
    α[i], V[i] = QGMM(CellDat[i], ones(size(CellDat[i])),
                        0.75, intercept = false)
end;
nCells = size(uqX, 1)
