%{
    This file performs all analayses requested in part 1 of problem set 3
    for Econ 761

    Date created:  23 Oct 2021
    Last modified: 23 Oct 2021
    Author: Danny Edgel
%}
clear; clc

%% Load the data and functions
load('data/ps2.mat'); load('data/iv.mat')

addpath('downloaded/code')
addpath('functions')

%% Estimate specifications

% OLS
est_ols = ols(s_jt, x1(:, 1), 1);

% OLS with brand FE
est_fe  = ols(s_jt, x1, 1);

% IV without brand FE
est_iv = ivreg(s_jt, x1(:, 1), iv(:, 2:end), 1);

% IV with brand FE
est_iv_fe = ivreg(s_jt, x1, [iv(:, 2:end), x1(:, 2:end)], 1);


% Print results
table = ['\\begin{tabular}{r|cccc}\n', ...
    '& OLS & OLS & IV & IV \\\\\\hline', ...
    '\\alpha & %4.3f & %4.3f  & %4.3f  & %4.3f \\\\', ...
    '& (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f) \\\\ &&&& \\\\', ...
    'FE? & & X & & X \\\\\\hline', ...
    '\\end{tabular}' ...
    ];

file = fopen('table1.tex', 'w');
fprintf(file, table, ...
    est_ols.b(1), est_fe.b(1), est_iv.b(1), est_iv_fe.b(1), ...
    est_ols.seb(1), est_fe.seb(1), est_iv.seb(1), est_iv_fe.seb(1) ...
    );
fclose(file);


%% Compute markups



%% Simulate post-merger equilibrium



%% Estimate mixed logit model

