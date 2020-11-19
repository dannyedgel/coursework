%%%
%%% This file is used to generate charts for Econ 712 problem set %2 of Q2
%%%
%%% Econ 712: Macroeconomics I
%%% Fall 2020
%%%
%%% Date created:  17 Nov 2020
%%% Last modified: 18 Nov 2020
%%% Author: Danny Edgel
%%%

%%% clear workspace
clear; clc

%%% set working directory
cd 'C:/Users/edgel/Google Drive/UW-Madison/f20/econ712/problem_sets/Q2PS2/'


%%%
%%% Question 3
%%%____________________________________________

%%% set VF convergence tolerance
tol = 1e-10;

%%% vector for number of scenarios (automating a, b, and c)
runs = [1,2,3]; %,2,3
gtitle = {'Policy function: Consumption', ...
    'Policy function: Consumption, \gamma = 1.01', ...
    'Transition dynamics: 20% shock to z'};

%%% define constants
beta  = 0.95;
delta = 0.1;
z     = 1;
gamma = 2;
alpha = 0.35; % Cobb-Douglas exponent on capital

%%% guess maximum capital level; define grid
maxkgrid = 5;
inkgrid = .01; % number of capital steps
minkgrid = 0;
kgrid  = minkgrid:inkgrid:maxkgrid;

n = length(kgrid); % number of elements in capital grid


%%% loop through scenarios
for run = runs
    %% print current run
    fprintf('\nRun %d\n',run);
    
    %% make current run's adjustment (defaults to 3(a))
    if run == 2; gamma = 1.01; old_ss = [k_ss,c_ss]; old_saddle = cpol; end
    if run == 3; gamma = 2; z = 1.2; end
   
    %% initialize loop variables
    v       = zeros(n);
    decis   = zeros(n);
    test    = 10; err = 1;
    iter    = 0;
    
    tic
    while err > tol || test ~= 0
        
        % calculate  value at each level of current capital, for each level
        % of next period capital
        cons = kgrid.^alpha + (1-delta)*kgrid - kgrid'*ones(1,n);
        util = (cons.^(1-gamma))/(1-gamma);
    
        % give large negative values to infeasible consumption choices
        % (i.e. enforce non-negativity constraint)
        util(cons<=0) = -10000;

        % sum the value function at the current capital level
        futval = sum(v,2);
        
        vint = util + beta*futval*ones(1,n);
        [t1,t2] = max(vint);
        
        tv     = t1; % RHS of Bellman
        tdecis = t2; % index of capital level that gave tv
        
        % calculate error and set up next iteration
        tv1     = tv';
        tdecis1 = tdecis';
        test= max(any(tdecis1-decis));
        err = max(max(abs(tv1-v)));
        nprin=10; % print every 10 iterations
        if iter==nprin*round(iter/nprin)
            disp(['Iteration = ', num2str(iter),... 
                ', k Error = ', num2str(max(max(abs(tdecis1-decis)))),...
                ', V Error = ', num2str(max(max(abs(tv1-v))))])
        end
        iter  = iter + 1;
        decis = tdecis1;
        v     = tv1;
        
    end
    toc
    
    %% calculate saddle path based on policy function
    kpol  = kgrid(decis);   % policy function for capital
    kpol  = kpol(2:end);
    kap   = kgrid(2:end);
    cpol  = z*kap.^alpha + (1-delta)*kap - kpol; % policy function for consumption
    
    %% define policy function for next period's capital in the consumption
    %% steady-state
    knext = ((1/alpha)*((1/beta)+delta-1))^(1/(alpha-1));
    
    %% define steady-state capital and consumption equations
    %% (e.g. delta c = 0)
    delk0 = z*kap.^alpha - delta*kap;
    delc0 = z*kap.^alpha + (1-delta)*kap - knext;
    k_ss = knext;
    c_ss = z*k_ss^alpha - delta*k_ss;
    
    %% plot phase diagram
    figure(run)
          plot(kap,delk0,kap,delc0,kap,cpol,'k--')
          hold on
          plot(k_ss,c_ss,'ko')
          if run > 1
              plot(kap,old_saddle,'r--')
              plot(old_ss(1),old_ss(2),'ro')
          end
          if run == 3
              quiver(old_ss(1),old_ss(2),k_ss-old_ss(1),c_ss-old_ss(2),.9);
          end
          hold off
          title(gtitle{run}); xlabel('K_t'); ylabel('C_t')
          if run == 1
              legend('\Delta k=0','\Delta c=0','Saddle Path', ...
                  'Steady State', 'location','Southeast')
          end
          if run == 2
              legend('\Delta k=0','\Delta c=0','Saddle Path', ...
                  'Steady State','Old Saddle Path', ...
                  'Old Steady State','location','Southeast')
          end
          if run == 3
              legend('\Delta k=0','\Delta c=0','Saddle Path', ...
                  'Steady State','Old Saddle Path', ...
                  'Old Steady State', 'One-time jump', ...
                  'location','Southeast')
          end
          saveas(gcf,['figure',num2str(run),'.png'])

end %% end current run
