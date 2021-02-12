%{
    This file is used to conduct all necessary tasks for problem set 3 of
    the first quarter of Econ 714

    Date created:  09 Feb 2021
    Last modified: 11 Feb 2021
    Author: Danny Edgel
%}

% clean workspace
clc; clear

% set convergence variables
tol             = 1e-5;  % convergence tolerance

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
    v.K_SS.^(-alpha) .* v.A_SS;

v.tauLt = ((1-v.TauL_SS)./v.TauL_SS).*( v.at + alpha*v.kt ...
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

% output LaTeX table of persistence parameters
filename = 'table4.tex';

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

tbl = ...
    ['\\begin{tabular}{r c}\n',...
    '\\hline $\\rho_{a}$ & %4.3f \\\\ \n',...
    '$\\rho_{g}$ & %4.3f \\\\ \n', ...
    '$\\rho_{\\tau_L}$ & %4.3f \\\\ \\hline\n',...
    '\\end{tabular}'];
fprintf(file1,tbl,...
    round(v.at_rho,3), round(v.gt_rho,3),round(v.tauLt_rho,3));
fclose(file1);


%% Question 5: Implement Blanchard-Kahn to solve the model

% assume persistence parameter for tau_It = 0
v.tauIt_rho = 0;

% calculate steady state of TauIt
v.TauI_SS = ((alpha*beta*v.A_SS.*(v.K_SS.^(alpha-1)).*(v.L_SS.^(1-alpha))) ...
    /(1-beta*(1-delta))) - 1;

% calculate steady state of X, which is equal to everything in the
% parentheses of the RHS of the Euler equation
X_SS = alpha*v.A_SS.*(v.K_SS.^(alpha-1)).*(v.L_SS.^(1-alpha)) ...
    + (1-delta)*(1 + v.TauI_SS);

% define theta and gamma to simplify later matrix construction
theta = alpha*v.A_SS.*(v.K_SS.^(alpha-1)).*(v.L_SS.^(1-alpha));
gamma = (theta*(alpha-1))./(X_SS*(phi+alpha));

% define A matrix
A = [ mean((alpha*delta*v.Y_SS)./v.I_SS + 1-delta + ...
    (delta*(1-alpha)*alpha*v.Y_SS)./(v.I_SS*(phi+alpha))), ...
    mean((delta*(alpha-1)*sigma*v.Y_SS)./(v.I_SS*(phi+alpha))-...
    (delta*v.C_SS)./v.I_SS); ...
    mean(((phi*gamma)./(1-gamma)).*( (alpha*delta*v.Y_SS)./(v.I_SS) ...
    +1-delta+ (delta*(1-alpha)*alpha*v.Y_SS)./(v.I_SS*(phi+alpha)) )), ...
    mean( (1./(1-gamma)).*( sigma + phi*gamma.* ...
    ( (delta*(alpha-1)*sigma*v.Y_SS)./(v.I_SS*(phi+alpha)) - ...
    (delta*v.C_SS)./(v.I_SS) ) ) ) ];

% return Q and Lambda matrices 
[Q, Lambda] = eig(A);
Qinv        = inv(Q);

% use explosive eigenvalue to calculate saddle path
if abs(Lambda(1,1)) > 1; r = 1; end % r gives the row of Qinv that 
if abs(Lambda(2,2)) > 1; r = 2; end % corresponds to the explosive lambda


c_saddle    = -(Qinv(r,1)/Qinv(r,2))*v.kt;

figure(1)
    plot(v.kt,c_saddle); yline(0); xline(0)
    xlabel('k_t'); ylabel('c_t'); axis equal
    title('Blanchard-Kahn Saddle Path')
    saveas(gcf,'figure5.png')
    
% save saddle ratio separately, for use in next problem
Q_rat = -(Qinv(r,1)/Qinv(r,2));

% output the mean, min, and max values of lower-left entry of A
a11 = (alpha*delta*v.Y_SS)./v.I_SS + 1-delta + ...
    (delta*(1-alpha)*alpha*v.Y_SS)./(v.I_SS*(phi+alpha));
a12 = (delta*(alpha-1)*sigma*v.Y_SS)./(v.I_SS*(phi+alpha))-...
    (delta*v.C_SS)./v.I_SS;
a21 = ((phi*gamma)./(1-gamma)).*( (alpha*delta*v.Y_SS)./(v.I_SS) ...
    +1-delta+ (delta*(1-alpha)*alpha*v.Y_SS)./(v.I_SS*(phi+alpha)) );
a22 = (1./(1-gamma)).*( sigma + phi*gamma.* ...
    ( (delta*(alpha-1)*sigma*v.Y_SS)./(v.I_SS*(phi+alpha)) - ...
    (delta*v.C_SS)./(v.I_SS) ) ) ;

filename = 'q5_a21.tex';

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

fprintf(file1, ...
    ['The mean of $a_{11}$, for example, is %4.5f, with a minimum and', ...
    ' maximum of %4.5f and %4.5f, respectively.'], ...
    round(mean(a11),5), round(min(a11),5), round(max(a11),5));
fclose(file1);


%% Question 6: Solve the fixed-point problem for tauIt

% initial conjecture for tauIt_rho: .8 (no reason for this value)
v.tauIt_rho = 0.8;

% in each period, calculate "cnext":= expected consumption in the next
% period, using the saddle path value of c, given this period's k_{t+1}
cnext = Q_rat*lagmatrix(v.kt,-1);

% calculate "lnext" as the expected labor in the next period, using
% expected consumption, etc. from the labor optimality condition
lnext = (v.at_rho + alpha*lagmatrix(v.kt,-1) - sigma*cnext ...
    - v.tauLt_rho*(v.TauL_SS./(1-v.TauL_SS)).*v.tauLt)/(phi + alpha);

% iteratively solve the model, estimate tauIt_rho, and repeat until the
% estimated value of tauIt_rho converges

diff    = 1;    % initialize convergence gap
attempt = 1;    % loop counter
while diff > tol && attempt < 200
    
    % infer tauIt from Euler eqn using current loop's value of tauIt_rho
    v.tauIt = ( sigma*X_SS.*v.TauI_SS.*(v.ct - cnext) + ...
        alpha*v.TauI_SS.*v.A_SS.*(v.K_SS.^(alpha-1)).*(v.L_SS.^(1-alpha)) ...
        .* (v.at_rho*v.at + (alpha-1)*lagmatrix(v.kt,-1) ...
        + (1-alpha)*lnext) ) ./(X_SS - v.TauI_SS*v.tauIt_rho*(1-delta));
    
    % re-estimate tauIt_rho and update the convergence gap
    tauIt_rho_new = regress(v.tauIt, lagmatrix(v.tauIt, 1));
    
    diff        = abs(tauIt_rho_new - v.tauIt_rho);
    v.tauIt_rho = tauIt_rho_new;
    
    fprintf('Covergence attempt %d; diff = %d\n', attempt, diff)
    attempt     = attempt + 1;
    
end

% output LaTeX code for persistence parameter
filename = 'q6.tex';

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

fprintf(file1,'${\\rho_{\\tau_I}=%4.3f}$', round(v.tauIt_rho,3));
fclose(file1);

%% Question 7: Plot dynamics of wedges
figure(2)
    plot(date1, v.at, date1, v.gt, date1, v.tauLt, date1, v.tauIt)
    datetick('x','yyyy'); yline(0)
    xlim([datenum('01/01/1980') datenum('12/31/2020')])
    title('Wedge Dynamics: deviations from steady-state')
    set(legend('$a_t$', '$g_t$', '$\hat{\tau}_{Lt}$','$\hat{\tau}_{It}$'),...
        'Interpreter', 'latex')
    saveas(gcf,'figure7.png')


%% Question 8: Solve for each wedge separately; plot output counterfactual

% define kappa separately to ease matrix construction
kappa = (delta*(1-alpha)*v.Y_SS)./((phi+alpha)*v.I_SS);

% Define B matrix
B = [ mean(kappa), -mean((delta*v.Y_SS)./(3*v.I_SS)), 0, ...
    -mean((v.TauL_SS./(1-v.TauL_SS)).*kappa); ...
    mean( (((1-alpha)*theta)./(sigma*(phi+alpha)*X_SS)).* ...
    (v.at_rho - (delta*(1-alpha)*v.Y_SS)./(v.I_SS) ) ), ...
    -mean( ((alpha-1)*delta*theta.*v.Y_SS)./(3*sigma*X_SS.*v.Y_SS) ), ...
    mean( ( (((1-delta)*v.tauIt_rho)./X_SS) - 1./v.TauI_SS ).* ...
    ((1 + v.TauI_SS)/sigma) ), ...
    mean( (((1-alpha)*theta.*v.TauL_SS)./ ...
    (sigma*(phi+alpha)*X_SS.*(1-v.TauL_SS))).* ...
    (v.tauLt_rho + (delta*(alpha-1)*v.Y_SS)./(v.I_SS) ) ) ];

% Set up counterfactual: define a new time series for each wedge, assuming
% all other wedges = 0

wedges = {'at', 'gt','tauIt', 'tauLt'};

% solve for each wedge separately
Theta = (-1/Lambda(r,r))*B(r,:)/(eye(4) - (1/Lambda(r,r))*diag([v.at_rho, v.gt_rho, ...
    v.tauIt_rho, v.tauLt_rho]));

for i = 1:length(wedges)
    
    % save current wedge in 'w'
    w = wedges{i};
    
    % initialize capital, consumption, and labor series
    cf.(w).kt = zeros(length(v.ct),1); cf.(w).ct = cf.(w).kt; 
    cf.(w).lt = cf.(w).kt; 
    cf.(w).kt(1) = v.kt(1);
    
    % initialize wedges 
    for j = 1:length(wedges)
        if i == j; mult = 1; else; mult = 0; end
        cf.(w).(wedges{j}) = mult*v.(wedges{j});
    end
    z = [ cf.(w).(wedges{1}), cf.(w).(wedges{2}), ...
        cf.(w).(wedges{3}), cf.(w).(wedges{4}) ]; 
    
    % loop though periods, using the B-K solution to solve for ct, then the
    % log-linearized model equations to define forward every other variable
    for t = 1:length(v.ct)
        
        % consumption
        cf.(w).ct(t) = (1/Qinv(r,2))*( -1*Qinv(r,1)*cf.(w).kt(t) ...
            + Theta*z(t,:)');
        
        % labor 
        cf.(w).lt(t) = (1/(phi+alpha))*( cf.(w).at(t) + ...
            alpha*cf.(w).kt(t) - sigma*cf.(w).ct(t) - ...
            (v.TauL_SS(t)/(1-v.TauL_SS(t)))*cf.(w).tauLt(t) );
        
        % next period's capital
        if t < length(v.ct)
            cf.(w).kt(t+1) = (1-delta)*cf.(w).kt(t) + ...
                (v.Y_SS(t)/v.I_SS(t))*( cf.(w).at(t) + alpha*cf.(w).kt(t) ...
                + (1-alpha)*cf.(w).lt(t) - ...
                (v.C_SS(t)/v.Y_SS(t))*cf.(w).ct(t) - (1/3)*cf.(w).gt(t));
        end
        
    end % end t loop
    
    % save output counterfactual
    cf.(w).yt = cf.(w).at + alpha*cf.(w).kt + (1-alpha)*cf.(w).lt;

end % end wedge loop

% plot wedge counterfactuals 
figure(3)
    plot(date1, v.yt, date1, cf.at.yt, date1, cf.gt.yt, date1, cf.tauLt.yt, ...
        date1, cf.tauIt.yt)
    datetick('x','yyyy'); yline(0)
    xlim([datenum('01/01/1980') datenum('12/31/2020')])
    title('Wedge Dynamics: GDP counterfactuals')
    set(legend('Actual','$a_t$', '$g_t$', '$\hat{\tau}_{Lt}$', ...
        '$\hat{\tau}_{It}$'), 'Interpreter', 'latex')
    saveas(gcf,'figure8.png')




