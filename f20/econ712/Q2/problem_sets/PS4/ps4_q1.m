%{
                    Question 1
   -------------------------------------------------
%}
%%% (a) solve for the value function


%% define transition matrix
Q = [0.85, 0.15; 0.05, 0.95];

%% define household parameters
beta    = 0.95;
gamma   = 3;
N = 2;          % # of wage states
s=[.7,1.1];       % employment states

if q1c == 1
    phi  = 2;  % borrowing limit  
else
    phi  = 0;  % borrowing limit 
end

%% define firm parameters
alpha = 0.36;
delta = 0.08;

%% save parameters in struct for passing on to vf_iter()
params          = struct();
params.Q        = Q;
params.beta     = beta;
params.gamma    = gamma;
params.s        = s;
params.phi      = phi;
params.alpha    = alpha;
params.delta    = delta;
                                                              
%% form capital grid
maxkap = 15;                     % maximum value of capital grid  
minkap = -phi;                   % borrowing constraint
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points



%% Calculate equilibrium--simulate market


% guess capital demand in the steady-state
k_demand = 5;
    
k_diff     = 1; % difference between current and next period capital 
k_iter     = 0; % current iteration
max_iter = 10000; % stop after 10,000 iterations if no convergence
kgrid = kap;
kap = kap(3:end);
count = 0;
% attempt market equilibrium
while abs(k_diff) > tol && k_iter < max_iter && count < 3
    k_iter = k_iter + 1;
    disp(['Attempt ',num2str(k_iter),' at equilibrium'])
    
    % iterate value function to retrieve capital accumulation policy function
    [decis,dind] = vf_iter(params,kgrid,k_demand,tol);
    decis_ind = dind(3:end);
    
    % attempt convergence between capital supply and demand
    disp('  Calculating Ks...')
    
    % initialize arbitrary asset distribution
    adist0    = ones(nkap,N);
    adist0    = adist0 ./(sum(adist0(:)));
    err = 1; % initialize tolerance
    i   = 0;
    while (err > tol && i < 1000)
        adist = 0*adist0;

        % loop through labor states and starting asset amounts
        for n = 1:N
           for a = 1:nkap
               targ = dind(a,n);
               adist(dind(a,n),1) = adist(targ,1)+adist0(a,n)*Q(n,1);
               adist(dind(a,n),2) = adist(targ,2)+adist0(a,n)*Q(n,2);
           end
        end
        i = i+1;
        err = sum(abs((adist0(:)-adist(:))));
        adist0 = adist;
        fprintf('     Iteration %d, err = %d\n',i,err)
    end
    amarg = sum(adist,2); % marginal distribution of assets in each l state
    
    % calculate capital supply using the marginal distributions
    k_supply = sum(amarg.*kgrid');
    
    % count number of repeated differences
    if k_demand - k_supply == k_diff
        count = count + 1;
    else
        count = 0;
    end
    
    k_diff = k_demand - k_supply;
    %k_demand = (k_supply+k_demand)/2;
    %k_demand = k_supply;
    [kmin,kind] = min(abs(0.95*k_demand + 0.05*k_supply - kgrid));
    k_demand = kgrid(kind);
    
    % report on k_ss convergence every iteration
    disp(['k_ss attempt = ', num2str(k_iter),... 
            ', k_ss Error = ', num2str(max(max(k_diff))),...
            ', k_demand = ', num2str(k_demand)])

end

    
%% save and/or calculate steady-state values
k_ss = k_demand;
r = alpha*k_ss^(alpha-1);
w = ((k_ss^alpha)-r*k_ss);

% calculate implicit consumption policy function
[xs,ykap] = meshgrid(s,kap);
condecis = w.*xs+(1+r-delta).*ykap-decis;

%% output results to .tex file
if q1c == 0
    filename = '1a.tex';
else
    filename = '1c.tex';
end

if exist(filename, 'file')==2
  delete(filename);
end
file1 = fopen(filename,'w');

tbl = ...
    ['\\begin{center}\n\\begin{tabular}{r c}\n',...
    '\\hline Capital  & %4.3f \\\\ \n',...
    'Interest rate & %4.3f \\\\ \n', ...
    'Wage & %4.3f \\\\ \\hline\n',...
    '\\end{tabular}\n\\end{center}'];
fprintf(file1,tbl,...
    round(k_ss,3), round(r,3),round(w,3));
fclose(file1);

if q1c == 0
    %%% (b) plot equilibrium distributions of income and consumption

    %% calculate distributions using the distribution of capital

    % consumption
    c       = [condecis(:,1);condecis(:,2)];
    prob_c  = [adist(3:end,1);adist(3:end,2)];
    sortmat = sort([c,prob_c],1);
    c = sortmat(:,1); prob_c = sortmat(:,2);
    if length(unique(c))<length(c)
       uni_c = unique(c);
       uniprob_c = zeros(length(uni_c));
       for i = 1:length(uni_c)
           uniprob_c(i) = sum(prob_c(c == uni_c(i)));
       end
       c = uni_c; prob_c = uniprob_c;
    end

    % income
    inc      = w.*xs+(1+r-delta).*ykap;
    inc      = [inc(:,1);inc(:,2)];
    prob_inc = [adist(3:end,1);adist(3:end,2)];
    sortmat = sort([inc,prob_inc],1);
    inc = sortmat(:,1); prob_inc = sortmat(:,2);
    if length(unique(inc))<length(inc)
       uni_inc = unique(inc);
       uniprob_inc = zeros(length(uni_inc));
       for i = 1:length(uni_inc)
           uniprob_inc(i) = sum(prob_inc(inc == uni_inc(i)));
       end
       inc = uni_inc; prob_inc = uniprob_inc;
    end

    %% plot PMFs
    figure(1)
        bar(c,prob_c, 'BarWidth', 1)
        title('PMF of equilibrium consumption')
        xlabel('c_t'); ylabel('Mass')
        saveas(gcf,'figure1bi.png')

    figure(2)
        bar(inc,prob_inc, 'BarWidth', 1)
        title('PMF of equilibrium income')
        xlabel('Income at t'); ylabel('Mass')
        saveas(gcf,'figure1bii.png')
end
    


