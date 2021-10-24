%{
    This file performs all analayses requested in part 1 of problem set 3
    for Econ 761

    Date created:  23 Oct 2021
    Last modified: 24 Oct 2021
    Author: Danny Edgel
%}
clear; clc

%% Load the data and functions
load('data/ps2.mat'); load('data/iv.mat')

addpath('downloaded/code')
addpath('functions')


%% Prepare data for MNL estimation

% normalize by first brand
norm_id = floor(id(1)/100000);
new_id  = id((floor(id/100000)~=norm_id));
cq      = (id - (floor(id/100000)*100000));
new_cq  = (new_id - (floor(new_id/100000)*100000));
X = full([x1(floor(id/100000) ~= norm_id, 1), ...
    x1(floor(id/100000) ~= norm_id, 3:end)]);
Z = iv(floor(id/100000) ~= norm_id, 2:end);
S = zeros(size(X, 1), 1);

for i = 1:length(unique(cq))
    new_i   = new_cq == cq(i);
    old_i   = cq == cq(i) & floor(id/100000) ~= norm_id;
    norm_i  = cq == cq(i) & floor(id/100000) == norm_id;
    X(new_i, 1) = X(new_i, 1) - x1(norm_i, 1);
    S(new_i, 1) = log(s_jt(old_i)) - log(s_jt(norm_i));
end


%% Estimate specifications

% OLS
est_ols = ols(S, X(:, 1), 1);

% OLS with brand FE
est_fe  = ols(S, X, 1);

% IV without brand FE
est_iv = regressIV(S, X(:, 1), [], Z, 1);

% IV with brand FE
est_iv_fe = regressIV(S, X(:, 1), X(:, 2:end), Z, 1);


% Print results
table = ['\\begin{tabular}{r|cccc}\n', ...
    '& OLS & OLS & IV & IV \\\\\\hline \n', ...
    '\\alpha & %4.3f & %4.3f  & %4.3f  & %4.3f \\\\ \n', ...
    '& (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f) \\\\ \n &&&& \\\\ \n', ...
    'FE? & & X & & X \\\\ \n &&&& \\\\ \n', ... 
    '$R^2$ & %3.2f & %3.2f & %3.2f & %3.2f \\\\\n', ...
    'N & %3.0f & %3.0f & %3.0f & %3.0f \\\\\\hline \n', ...
    '\\end{tabular}' ...
    ];

file = fopen('table1.tex', 'w');
fprintf(file, table, ...
    est_ols.b(1), est_fe.b(1), est_iv.b(1), est_iv_fe.b(1), ...
    est_ols.seb(1), est_fe.seb(1), est_iv.seb(1), est_iv_fe.seb(1), ...
    est_ols.R2(1), est_fe.R2(1), est_iv.R2(1), est_iv_fe.R2(1), ...
    size(S, 1), size(S, 1), size(S, 1), size(S, 1) ...
    );
fclose(file);


%% Compute markups



%% Simulate post-merger equilibrium



%% Estimate mixed logit model

