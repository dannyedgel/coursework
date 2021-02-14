%{
    This file is used to perform the tasks required for question 3 of
    problem set 3 for the first quarter of ECON 710

    Date created:  13 Feb 2021
    Last modified: 13 Feb 2021
    Author: Danny Edgel
%}

% clean workspace
clc; clear

%% Step 1: Import and prepare Angrist & Krueger (1991) data

% read data from .csv file
dat = readtable('AK91.csv');

% Extract X1 and vectors
Y  = dat.lwage;
X1 = dat.educ;

% generate X2 matrix -- for state and year of birth dummies, remove states 
% and years with no observations (otherwise X won't be invertible)
v = {'sob', 'yob'};

for i = 1:length(v)
    x.(v{i})    = dummyvar(dat.(v{i}));
    count       = sum(x.(v{i}));
    x.(v{i})    = x.(v{i})(:, count>0);
end

% (exclude first dummy of each dummy matrix)
X2 = [ ones(size(X1)), x.yob(:,1:(end-1)), x.sob(:,1:(end-1)) ]; 

% save X separately, for coding ease in step 2
X = [ X1, X2 ];

% generate Z matrix from qob and X2 (exclude first quarter)
qob_dummies = dummyvar(dat.qob);

Z = [ qob_dummies(:,2:4), X2 ];

% save sample size
n = size(X1, 1);

%% Step 2: calculate 2SLS estimator for beta_1 and its heteroskedasticity-
%   robust standard error from (12.40) in Hansen's textbook

% estimate beta
bhat =  (X'*Z/(Z'*Z)*Z'*X)\(X'*Z/(Z'*Z)*Z'*Y);

% estimate the constituent matrices of the standard error
Q_zz = (1/n)*(Z'*Z);
Q_xz = (1/n)*(X'*Z);
ehat = Y - X*bhat;
Ohat = 0*Q_zz;
for i=1:n; Ohat = Ohat + (1/n)*Z(i,:)'*Z(i,:)*ehat(i)^2; end

% calculate standard error
Vb = ((Q_xz/Q_zz*Q_xz')\(Q_xz/Q_zz*Ohat/Q_zz*Q_xz')/(Q_xz/Q_zz*Q_xz'))/n;


%% Step 3: output results

% initialize LaTeX file
filename = 'q3.tex';

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

% output bhat_1 and its heteroskedasticity-robust SE
tex = ...
    ['\\begin{align*}\n',...
    '\\beta^{2SLS}_1 &= %4.3f \\\\ \n',...
    '\\text{SE}      &= %4.3f  \n', ...
    '\\end{align*}'];
fprintf(file1,tex,...
    round(bhat(1), 3), round(sqrt(Vb(1, 1)), 3));
fclose(file1);