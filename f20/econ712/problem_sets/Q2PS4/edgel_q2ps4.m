%{
    This file is used to conduct all necessary tasks for problem set 4 of
    the second quarter of Econ 712

    Date created:  08 Dec 2020
    Last modified: 10 Dec 2020
    Author: Danny Edgel
%}

%%% clean workspace
clc; clear

%%% set VF iteration variables
tol = 1e-3;
inckap = 0.01;           % size of capital grid increments

%%% choose which files to run
q1 = 0;
q2 = 0;
q1c = 0;

%%% run files according to selection
if q1 == 1; ps4_q1; end
if q2 == 1; ps4_q2; end

%%% read table from 1a and 1c to create a comparative table

%% save file names and variable names (from column of LaTeX tables) in
%% cells for extracting values in loop
files   = {'1a.tex','1c.tex'};
vars    = {'\\hline Capital','Interest rate','Wage'};

%% loop through each file and variable to extract the numbers to be output
%% in the new table

values = cell(length(vars),length(files));  % contains properly-indexed
                                            % table values

for i = 1:length(files)
    % read existing table as single string
    atbl = fileread(files{i});
    
    for j = 1:length(vars)
        % extract line that contains current loop's variable
        vline = regexp(atbl,['[^\n]*',vars{j},'[^\n]*'],'match');
        
        % extract the variable's value from vline
        vstart      = regexp(vline{1},'&','end')+1;
        vend        = regexp(vline{1},'\\\\','start')-1;
        values{i,j} = vline{1}(vstart:vend);
    end
    
end



%% write new file comparing the two tables
filew = fopen('1_compare.tex','w');

tbl = ['\\begin{center}\n \\begin{tabular}{r c c}\n',...
    ' & $a_t\\geq0$ & $a_t\\geq-2$ \\\\ \\hline ',...
    'Capital & ',values{1,1},' & ',values{2,1},' \\\\',...
    'Interest rate & ',values{1,2},' & ',values{1,2},' \\\\',...
    'Wage & ',values{1,3},' & ',values{2,3},'\\\\ \\hline',...
    '\\end{tabular}\n \\end{center}'];

fprintf(filew,tbl);
fclose(filew);