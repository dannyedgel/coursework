% This function prints table 5 to table5.tex

function [] = printTable5(mu, p)


table = ['\\begin{tabular}{r|ccc}\n', ...
    ' & Mean & Median & Variance \\\\\\hline &&& \\\\ \n', ...
    ' Markups   & %4.3f & %4.3f & %4.3f \\\\ \n', ...
    ' MC        & %4.3f & %4.3f & %4.3f \\\\ \n', ...
    ' &&& \\\\\\hline \n', ...
    '\\end{tabular}' ...
    ];

outvec = [mean(mu); median(mu); var(mu); ...
            mean(p - mu); median(p-mu); var(p-mu)];

file = fopen('table5.tex', 'w');
fprintf(file, table, outvec);

fclose(file);
end