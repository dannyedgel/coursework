%{
    This file is used to conduct all necessary tasks for problem set 4 of
    the first quarter of Econ 714

    Date created:  20 Feb 2021
    Last modified: 20 Feb 2021
    Author: Danny Edgel
%}

%% Preamble: set up model primitives

% clean workspace
clc; clear

% set convergence variables
tol             = 1e-5; % convergence tolerance
max_j           = 1000;  % maximum convergence loops
delta           = .7;   % convergence tuning parameter

% set seed for reproducability
rng(4156);

% set an epsilon error for rho
err = 1e-8;

% set parameters
rho     = 1 + err;
theta   = 5;
Nk      = 20;
K       = 100000;
W       = 1; % (assumption)

% generate Aik matrix, with firms in rows and industries in columns
logAik  = normrnd(0, 1, Nk, K);
Aik     = exp(logAik);

%% Prepare model solution loop

% initialize s_ik with even distribution across all firms in all industries
s0  = (ones(Nk, K) / Nk).^(1-theta);

% initialize loop variables
j       = 0;
diff    = 1;

%% Solve model: simulate until fixed point is reached

while ( diff > tol && j < max_j )
    
    % update loop count
    j = j + 1;
    
    % calculate prices, given current loop's sik
    eta = 1 - (1./((theta-rho)*s0+1-theta)); % firm mark-up
    Pik = (W ./ Aik) .* eta;
    
    % calculate industry price indices
    Pk  = sum(Pik.^(1-theta), 1).^(1/(1-theta));
    
    % update s_ik
    sik = (Pik ./ Pk).^(1-theta);
    
    % calculate difference in s_ik values (update loop)
    diff    = norm(abs(sik-s0));
    s0      = delta*s0 + (1-delta)*sik;
    
    % print loop summary
    fprintf('\nLoop %d, difference: %12.5f', j, diff)
    
end
fprintf('\n')

% calculate real wage
P       = ((1/K)*sum(Pk.^(1-rho))).^(1/(1-rho));
w_real  = W / P;  

% calculate 'first-best' consumption
Pik_fb      = W./Aik;
Pk_fb       = sum(Pik_fb.^(1-theta), 1).^(1/(1-theta));
P_fb        = ((1/K)*sum(Pk_fb.^(1-rho))).^(1/(1-rho));
w_real_fb   = W / P_fb;

%% Output results for question 5

% initialize LaTeX file
filename = 'q5.tex';

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

% output real wages
fprintf(file1,'&C = %4.3f &C^{fb} = %4.3f',...
    round(w_real, 3), round(w_real_fb, 3));
fclose(file1);
