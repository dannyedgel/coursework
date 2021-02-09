%{
    This file is used to conduct all necessary tasks for problem set 3 of
    the first quarter of Econ 714

    Date created:  09 Feb 2021
    Last modified: 09 Feb 2021
    Author: Danny Edgel
%}

% clean workspace
clc; clear

% set convergence variables
tol             = 1e-3;  % convergence tolerance

%% Question 1: import and set up FRED data

% establish connection and set date ranges
url     = 'https://research.stlouisfed.org/fred2/';
con     = fred(url);
start1  = '01/01/1980';
start2  = '01/01/1950';
dateend = '12/31/2020';

Ct = fetch(con, 'PCECC96', start1, dateend).Data(:,2);
Lt = fetch(con, 'PAYEMS', start1, dateend).Data(:,2);
Yt = fetch(con, 'GDPC1', start1, dateend).Data(:,2);
It = fetch(con, 'GPDIC1', start2, dateend).Data(:,2);

% aggregate labor to quarterly
Lt_agg = zeros(length(Ct),1);
for t = 1:length(Ct)
    Lt_agg(t) = mean(Lt((3*t-2):(3*t)));
end
Lt = Lt_agg;

% save variable names and store in struct
varnames    = {'Ct', 'Lt', 'Yt', 'It'};
vars        = {Ct, Lt, Yt, It};

for i = 1:length(vars); v.(varnames{i}) = vars{i}; end

% save the index of It that corresponds to the first date of the other
% series
date1 = fetch(con, 'PCECC96', start1, dateend).Data(:,1);
date2 = fetch(con, 'GPDIC1', start2, dateend).Data(:,1);
date2 = [ [1:length(date2)]', date2 ];
q1_1980 = date2(date2(:,2) == date1(1), 1);

% test q1_1980 (need to do only once)
%datestr(date2(q1_1980,2))


%% Question 2: log and de-trend with Hodrick-Prescott

% loop through variables, taking logs and detrending
for i = 1:length(vars)
    
    % define the name of the current variable's steady state
    ss_name = [varnames{i}(1), '_SS'];
    
    
    % smooth using pre-set hpfilter function, storing variable in struct
    % (suppress first output, which is the trend; we want the deviation 
    % from the trend)
    [v.(ss_name), v.(lower(varnames{i}))] = ...
        hpfilter(log(vars{i}), 1600);
    
end



%% Question 3: Estimate capital stock

delta = 0.025;

% generate series of deviations from the capital steady state; assume 0 in
% 1950q1; after that, k_{t-1} = (1-delta)kt + it

v.kt = zeros(length(It),1);

for t = 2:length(It)
    v.kt(t) = (1-delta)*v.kt(t-1) + v.it(t-1);
end

% subset capital and investment vectors, saving long versions
v.kt_long = v.kt; v.kt = v.kt(q1_1980:end);
v.it_long = v.it; v.it = v.it(q1_1980:end);

v.I_SS_long = v.I_SS; v.I_SS = v.I_SS(q1_1980:end);

% calculate capital SS from law of motion
v.K_SS = (1/delta)*v.I_SS;



%% Question 4: Estimate additional model variables and compute persistence
%               parameters

% Define parameters
alpha   = 1/3;
sigma   = 1;
phi     = 1;
beta    = 0.99;

% estimate a_t from the linearized market clearing condition
v.at = (1/phi)*lagmatrix(v.kt,-1) - (1/(phi*beta))*v.kt ...
    + ((phi-delta)/phi)*v.ct;

v.A_SS = ones(length(v.at),1); % A_SS assumed = 1

% estimate g_t from log-linearized Yt = Ct + It + Gt
v.gt = 3*( v.at + alpha*v.kt + (1-alpha)*v.lt - ...
    (v.C_SS./v.Y_SS).*v.ct - (v.I_SS./v.Y_SS).*v.it );

% estimate tau_Lt from labor supply, beginning with its steady-state
v.TauL_SS = 1 + v.L_SS.^(phi+alpha) .* v.C_SS.^(-sigma) .* ...
    v.K_SS.^(-alpha) .* v.A_SS.^(-2);

v.tauLt = ((1-v.TauL_SS)./v.TauL_SS).*( 2*v.at + alpha*v.kt ...
    - (alpha+phi)*v.lt - sigma*v.ct );

% compute persistence parameter for each wedge
wedges = {'at', 'gt', 'tauLt'}; 

for i = 1:length(wedges)
    
    % save cureent wedge
    w = v.(wedges{i});
   
    % save name of persistence parameter for current wedge
    w_name = [wedges{i}, '_rho'];
    
    % obtain parameter from AR(1) regression
    v.(w_name) = regress(w, lagmatrix(w,1));
    
end


%% Question 5: Implement Blanchard-Kahn to solve the model

% assume persistence parameter for tau_It = 0
v.tauIt_rho = 0;

% define eta to simplify later matrix construction
eta = beta*(1-alpha)*alpha*phi*(phi-delta)*(sigma^(-1));

% define A and B matrices
A = [ 1/beta, -(phi-delta); -((1-alpha)*alpha*phi/sigma), 1 + eta];
B = [ phi; ((beta*alpha*phi)/sigma)*(v.tauLt_rho - (1-alpha)*phi) ];

% return Q and Lambda matrices 
[Q, Lambda] = eig(A);
Qinv        = inv(Q);

% use explosive eigenvalue to calculate saddle path
if abs(Lambda(1,1)) > 1; r = 1; end % r gives the row of Qinv that 
if abs(Lambda(2,2)) > 1; r = 2; end % corresponds to the explosive lambda

c_saddle = -(Qinv(r,1)/Qinv(r,2))*v.kt;

figure(1)
    plot(v.kt,c_saddle); yline(0); xline(0)
    xlabel('k_t'); ylabel('c_t'); axis equal
    title('Blanchard-Kahn Saddle Path')
    saveas(gcf,'figure5.png')


%% Question 6: Solve the fixed-point problem for tauIt

% initial conjecture for tauIt_rho: .8 (no reason for this value)
v.tauIt_rho = 0.8;

% iteratively solve the model, estimate tauIt_rho, and repeat until the
% estimated value of tauIt_rho converges

diff = 1;   % initialize convergence gap

while diff > tol
    
    % solve model for consumption as a function of capital and wedges
    
end

%% Question 7: Plot dynamics
%{
figure(1)
    plot(t_vec,k_path,t_vec,k_ss*ones(151),'r--')
    text(t_vec(end)*.8,k_ss*1.0003,'Steady State','Color','red')
    xline(12,'--'); text(15,(max(k_path)+min(k_path))/2,'Event')
    title('Capital post-shock transition path')
    xlabel('t'); ylabel('K_t')
    saveas(gcf,'figure4a.png')
    
figure(2)
    plot(t_vec,c_path,t_vec,c_ss*ones(151),'r--')
    text(t_vec(end)*.8,c_ss*.9999,'Steady State','Color','red')
    xline(12,'--'); text(15,(max(c_path)+min(c_path))/2,'Event')
    title('Consumption post-shock transition path')
    xlabel('t'); ylabel('C_t')
    saveas(gcf,'figure4b.png')
%}

%% Question 8: 










