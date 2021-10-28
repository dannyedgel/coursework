% This function prints table 2 to table2.tex

function [] = printTable2(res, names)
markups = res;
table = ['\\begin{tabular}{r|cccc}\n', ...
    ' & (1) & (2) & (3) & (4) \\\\\\hline &&&& \\\\ \n', ...
    '$\\E{\\mu_{jt}}$ & %4.3f & %4.3f  & %4.3f  & %4.3f \\\\ \n', ...
    'Med($\\mu_{jt}$)& %4.3f & %4.3f & %4.3f & %4.3f \\\\\n ', ...
    'Var($\\mu_{jt}$)& %4.3f & %4.3f & %4.3f & %4.3f \\\\\n &&&&\\\\ \n', ...
    '$\\E{c_{jt}}$ & %4.3f & %4.3f  & %4.3f  & %4.3f \\\\ \n', ...
    'Med($c_{jt}$)& %4.3f & %4.3f & %4.3f & %4.3f \\\\\n ', ...
    'Var($c_{jt}$)& %4.3f & %4.3f & %4.3f & %4.3f \\\\\n &&&&\\\\ \n', ...
    '$\\E{m_{jt}}$ & %4.3f & %4.3f  & %4.3f  & %4.3f \\\\ \n', ...
    'Med($m_{jt}$) & %4.3f & %4.3f & %4.3f & %4.3f \\\\\n', ...
    'Var($m_{jt}$) & %4.3f & %4.3f & %4.3f & %4.3f \\\\\n', ...
    '&&&&\\\\\\hline \n', ...
    '\\end{tabular}' ...
    ];

file = fopen('table2.tex', 'w');
fprintf(file, table, ...
    mean(markups.(names{1}).mu), mean(markups.(names{2}).mu), ...
        mean(markups.(names{3}).mu), mean(markups.(names{4}).mu),  ...
    median(markups.(names{1}).mu), median(markups.(names{2}).mu), ...
        median(markups.(names{3}).mu), median(markups.(names{4}).mu),  ...
    var(markups.(names{1}).mu), var(markups.(names{2}).mu), ...
        var(markups.(names{3}).mu), var(markups.(names{4}).mu), ...
    mean(markups.(names{1}).mc), mean(markups.(names{2}).mc), ...
        mean(markups.(names{3}).mc), mean(markups.(names{4}).mc),  ...
    median(markups.(names{1}).mc), median(markups.(names{2}).mc), ...
        median(markups.(names{3}).mc), median(markups.(names{4}).mc),  ...
    var(markups.(names{1}).mc), var(markups.(names{2}).mc), ...
        var(markups.(names{3}).mc), var(markups.(names{4}).mc), ...    
    mean(markups.(names{1}).margins), mean(markups.(names{2}).margins), ...
        mean(markups.(names{3}).margins), ...
        mean(markups.(names{4}).margins),  ...
    median(markups.(names{1}).margins), ...
        median(markups.(names{2}).margins), ...
        median(markups.(names{3}).margins), ...
        median(markups.(names{4}).margins),  ...
    var(markups.(names{1}).margins), var(markups.(names{2}).margins), ...
        var(markups.(names{3}).margins), ...
        var(markups.(names{4}).margins) ...
    );
fclose(file);

end