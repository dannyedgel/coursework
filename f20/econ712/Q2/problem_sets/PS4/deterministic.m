function [kpol, cpol, k_ss, c_ss] = deterministic(params,kgrid)

%%% set VF convergence tolerance
tol = 1e-10;

%%% define constants
beta  = params.beta;
delta = params.delta;
z     = params.z;
gamma = params.gamma;
alpha = params.alpha; % Cobb-Douglas exponent on capital

n = length(kgrid); % number of elements in capital grid


%% initialize loop variables
v       = zeros(n);
decis   = zeros(n);
test    = 10; err = 1;
iter    = 0;

tic
while err > tol || test ~= 0

    % calculate  value at each level of current capital, for each level
    % of next period capital
    cons = z*kgrid.^alpha + (1-delta)*kgrid - kgrid'*ones(1,n);
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
kpol  = kpol(3:end);
kap   = kgrid(3:end);
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

end