%{
    This file is used to conduct all necessary tasks for problem set 4 of
    the second quarter of Econ 712

    Date created:  08 Dec 2020
    Last modified: 08 Dec 2020
    Author: Danny Edgel
%}

%%% clean workspace
clc; clear

%%% set VF iteration variables
tol = 1e-3;
inckap = 0.01;           % size of capital grid increments

%%% choose which files to run
q1 = 1;
q2 = 0;
q1c = 0;

%%% run files according to selection
if q1 == 1; ps4_q1; end
if q2 == 1; ps4_q2; end

%%% read table from 1a and 1c to create a comparative table

%% read and parse table 1a

% read table 1a as a string of text
atbl = fileread('1a.tex');

% extract equilibrium values for each variable
tkap = regexp(atbl,'[^\n]*\\hline Capital[^\n]*','match');
tr = regexp(atbl,'[^\n]*Interest rate[^\n]*','match');
tw = regexp(atbl,'[^\n]*Wage[^\n]*','match');

tkap1 = regexp(tkap{1},'&','end')+1;
tkap2 = regexp(tkap{1},'\\\\','start')-1;
tkap = tkap{1}(tkap1:tkap2);

tr1 = regexp(tr{1},'&','end')+1;
tr2 = regexp(tr{1},'\\\\','start')-1;
tr = tr{1}(tr1:tr2);

tw1 = regexp(tw{1},'&','end')+1;
tw2 = regexp(tw{1},'\\\\','start')-1;
tw = tw{1}(tw1:tw2);

values_1a = {tkap,tr,tw};

%% read and parse table 1c

% read table 1a as a string of text
atbl = fileread('1c.tex');

% extract equilibrium values for each variable
tkap = regexp(atbl,'[^\n]*\\hline Capital[^\n]*','match');
tr = regexp(atbl,'[^\n]*Interest rate[^\n]*','match');
tw = regexp(atbl,'[^\n]*Wage[^\n]*','match');

tkap1 = regexp(tkap{1},'&','end')+1;
tkap2 = regexp(tkap{1},'\\\\','start')-1;
tkap = tkap{1}(tkap1:tkap2);

tr1 = regexp(tr{1},'&','end')+1;
tr2 = regexp(tr{1},'\\\\','start')-1;
tr = tr{1}(tr1:tr2);

tw1 = regexp(tw{1},'&','end')+1;
tw2 = regexp(tw{1},'\\\\','start')-1;
tw = tw{1}(tw1:tw2);

values_1c = {tkap,tr,tw};

%% write new file comparing the two tables
filew = fopen('1_compare.tex','w');

tbl = ['\\begin{center}\n \\begin{tabular}{r c c}\n',...
    ' & $a_t\\geq0$ & $a_t\\geq-2$ \\\\ \\hline ',...
    'Capital & ',values_1a{1},' & ',values_1c{1},' \\\\',...
    'Interest rate & ',values_1a{2},' & ',values_1c{2},' \\\\',...
    'Wage & ',values_1a{3},' & ',values_1c{3},'\\\\ \\hline',...
    '\\end{tabular}\n \\end{center}'];

fprintf(filew,tbl);
fclose(filew);