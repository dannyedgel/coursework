function [kpol, ind] = vf_iter(params,kap,k_demand, tol)
%
% This function performs value function iterations over a capital grid,
% given capital demand and a set of parameters and iteration rules, to
% output capital and consumption policy functions
%
% Inputs:
%   params: a struct of (named) model parameters
%   kap: a capital grid
%   k_demand: constant capital demand
%   tol: tolerance level for VF convergence
%   
% Outputs:
%   kpol: capital policy function matrix (# of labor states in columns)
%   ind: index of each kprime value in kpol on the capital grid
%
% Date created:  09 Dec 2020
% Last modified: 09 Dec 2020
% Author: Danny Edgel
% 

%% extract parameters from the params struct
Q       = params.Q;
beta    = params.beta;
gamma   = params.gamma;
s       = params.s;
phi     = params.phi;
alpha   = params.alpha;
delta   = params.delta;
N       = length(s);          % # of wage states 

%  initialize some variables
nkap    = length(kap);
v       = zeros(nkap,N);
decis   = zeros(nkap,N);

test    = 10;

%  iterate on Bellman equation and get the decision 
%  rules and the value function at the optimum         
%
cons = zeros(nkap,nkap,N);
util = zeros(nkap,nkap,N);
vint = zeros(nkap,nkap,N);
tv   = zeros(N,nkap);
tdecis = zeros(N,nkap);

iter=0;
err = 1;

tic
while  test ~= 0 % || err > tol


    %  tabulate the utility function such that for zero or negative
    %  consumption utility remains a large negative number so that
    %  such values will never be chosen as utility maximizing      
    for j=1:N
        
        % calculate wages and interest rates for each capital decision
        r = alpha*k_demand^(alpha-1);
        w = (k_demand^alpha)-r*k_demand;
        
        % calculate consumption and utility
        cons(:,:,j) = s(j).*w + (1+r-delta).*(ones(nkap,1)*kap) -...
            kap'*ones(1,nkap);
        
        util(:,:,j) = (cons(:,:,j).^(1-gamma))/(1-gamma);
        
        % assign large negative values to infeasible consumption choices
        utilj = util(:,:,j);
        
        Aj = cons(:,:,j);                % give negative value to 
        i = find( Aj <= 0);              % infeasible consumption choice
        utilj(i) = -Inf;
        util(:,:,j) = utilj;
        
        
        futval = sum((ones(nkap,1)*Q(j,:)).*v,2); % nkap by 1 vector of
                                                  % future value function
                                                  
        % calculate RHS of the Bellman for all possible next-period capital
        % choices and select the maximum
        vint(:,:,j) = util(:,:,j) + beta*futval*ones(1,nkap);
        vint(i) = -Inf;
        
        [t1,t2] = max(vint(:,:,j));
        tv(j,:)     = t1; % maximum value
        tdecis(j,:) = t2; % index of next-period capital that gave max val
        
    end % end labor states loop
    
    tv1     = tv';
    tdecis1 = tdecis';
    
    % update convergence metrics
    test=max(any(tdecis1-decis));
    diff = tv1-v;
    err = norm(diff(~isnan(diff)));
    
    % report on convergence every 10 iterations
    nprin=10;
    if iter==nprin*round(iter/nprin)
        disp(['     Iteration=', num2str(iter),... 
            ', k Error=', num2str(max(max(abs(tdecis1-decis)))),...
            ', V Error=', num2str(err)])
    end
    
    % prepare for next iteration
    iter    = iter+1;
    decis   = tdecis1;
    v       = tv1;
    
end % end outermost while loop
toc

% focus on decisions off lowest boundary...
ind     = decis;              % index policy function
kpol    = -phi + (decis(3:end,:)-1)*(kap(2)-kap(1));     % actual policy function


end