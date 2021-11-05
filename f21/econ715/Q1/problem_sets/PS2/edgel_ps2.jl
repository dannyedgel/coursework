#==
    This file completes all of the tasks for 1(e) of problem
    set 2 for Econ 715

    Date created:  04 November 2021
    Last modified: 05 November 2021
    Author: Danny Edgel
==#


using Printf

# load functions
include("./econ715/Q1/problem_sets/PS2/functions.jl")

# simulate data with n = 100 and s (number of simulations) = 1000
param, sim = Initialize(100; s = 1000);
@unpack S, s = sim

# Initialize results arrays
θhat  = zeros(size(S, 3));
se    = zeros(6, size(S, 3))

# save Latex label of each SE estimate
labels = ["\$s.e.^{asy}\$", "\$s.e.^{bt}\$", "\$s.e.^{IQR}\$",
            "\$s.e.^{bt-r}\$", "\$s.e.^{sbt-normal}\$", 
            "\$s.e.^{sbt-rad}\$"];

for i = 1:length(θhat)

    if i % 100 == 0
        println("Simulation $i of $s...\n")
    end

    # save X and Y 
    Y, X = S[:, 1, i], S[:, 2, i];

    # estimate θ for each sample
    θhat[i]            = OLS(Y, X);
    se[1, i]           = AsympSE(Y, X);
    se[2, i], se[3, i] = BootstrapSE(Y, X);
    se[4, i], _        = BootstrapSE(Y, X; modify = true);
    se[5, i]           = ScoreBootstrap(Y, X);
    se[6, i]           = ScoreBootstrap(Y, X; binomial = true);
end

# write latex table of results
fname = "./econ715/Q1/problem_sets/PS2/table1.tex";
open(fname, "w") do io
    write(io, "\\begin{tabular}{r|ccc}\n")
    write(io, " & Mean & SD & Coverage \\\\\\hline &&& \\\\ \n")
    
    for i = 1:6
        b, c, d = seStats(θhat[:], se[i, :], param)
        str = @sprintf " & %1.3f & %1.3f & %1.3f \\\\ \n" b c d
        write(io, labels[i]*str)
    end

    write(io, " &&& \\\\\\hline \n \\end{tabular}")
end;
