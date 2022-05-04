#==
    This file performs all necessary tasks for problem set 2
    of the second quarter of Econ 717

    Date created:  02 May 2022
    Last modified: 04 May 2022
    Author: Danny Edgel (edgel@wisc.edu)
==#


### Housekeeping

# load necessary packages
using Optim, Random, Plots, Printf, Combinatorics

# declare user-written functions 
include("./ps2b_q1functions.jl")
#include("./ps2b_q2functions.jl")


### Question 1

# initialize empty arrays for table 2
M, β, V, F = [], [], [], []

## a) Write a program that simulates data for 2,000 observations
dat = SimulateData(2000)

## b) obtain OLS point estimates for a model with sᵢ and an intercept 
res = OLS(dat[:, 1], dat[:, 2])

push!(M, "OLS")
push!(β, res[1])
push!(V, res[2])
push!(F, Inf)

## c) test each instrument, calculating the correlation coefficient
corr = cor(dat[:, 2:end])
str = @sprintf "\\begin{tabular}{rc}
    Instrument & Correlation Coefficient \\\\\\hline
    \$z\\_1\$ & %0.3f \\\\
    \$z\\_2\$ & %0.3f \\\\
    \$z\\_3\$ & %0.3f
    \\end{tabular}" corr[2, 1] corr[3, 1] corr[4, 1]

open("q1c.tex", "w") do io
    write(io, str)
end

## d) compute the 2SLS estimate for β using all possible combinations
##    of instruments, reporting the estimates, SEs, and first-stage
##    F-statistic
for x in collect(combinations(3:size(dat, 2)))
    res = TwoSLS(dat[:, 1], dat[:, 2], dat[:, x], Ftest=true)
    push!(M, join(x .- 2, ", "))
    push!(β, res[1])
    push!(V, res[2])
    push!(F, res[3])
end

## f) repeat a-d with 500,000 observations and print all results 
dat = SimulateData(500000)

res = OLS(dat[:, 1], dat[:, 2])

push!(M, "OLS")
push!(β, res[1])
push!(V, res[2])
push!(F, Inf)

for x in collect(combinations(3:size(dat, 2)))
    res = TwoSLS(dat[:, 1], dat[:, 2], dat[:, x], Ftest=true)
    push!(M, join(x .- 2, ", "))
    push!(β, res[1])
    push!(V, res[2])
    push!(F, res[3])
end


# print results to a LaTeX table
k = Int(length(β) / 2)
str = @sprintf "\\begin{tabular}{r|cccccc} 
    \\hline\\hline
    & \\multicolumn{3}{c}{\$N = 2,000\$} & \\multicolumn{3}{c}{\$N = 500,000\$} \\\\
    & \$\\beta_1\$ & \$\\beta_2\$ & F-Stat & \$\\beta_1\$ & \$\\beta_2\$ & F-Stat \\\\\\hline\\\\
    "

for i in 1:k
    str = join([str
        @sprintf "%s & %2.2f & %2.2f & %6.1f & " M[i] β[i][1] β[i][2] F[i]
    ])
    str = join([str
        @sprintf "%2.2f & %2.2f & %6.1f \\\\
        " β[i+k][1] β[i+k][2] F[i+k]
    ])
    str = join([str
        @sprintf " & (%2.2f) & (%2.2f) &  & " sqrt(V[i][1, 1]) sqrt(V[i][2, 2])
    ])
    str = join([str
        @sprintf "(%2.2f) & (%2.2f) &  \\\\ \\\\
        " sqrt(V[i+k][1, 1]) sqrt(V[i+k][2, 2])
    ])
end
str = join([str
    @sprintf "\\\\\\hline\\hline 
    \\end{tabular}"])

open("table2.tex", "w") do io
    write(io, str)
end
