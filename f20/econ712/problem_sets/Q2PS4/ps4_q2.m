%{
                    Question 2
   -------------------------------------------------
%}
%%% define constants
beta  = 0.95;
delta = 0.1;
gamma = 2;
alpha = 0.35; % Cobb-Douglas exponent on capital

%%% define stochastic variables
z = [0.8,1.2];          % productivity states
P = [0.9,0.1;0.1,0.9];  % transition matrix
N = length(z);

%%% guess maximum capital level; define grid
maxkgrid = 10;
inkgrid = inckap; % number of capital steps
minkgrid = 0;
kgrid  = minkgrid:inkgrid:maxkgrid;

n = length(kgrid); % number of elements in capital grid

%% initialize loop variables
v       = zeros(n,N);
decis   = zeros(n,N);
cons = zeros(n,n,N);
util = zeros(n,n,N);
vint = zeros(n,n,N);
tv   = zeros(N,n);
tdecis = zeros(N,n);

test    = 10; err = 1;
iter    = 0;
tic
while err > tol || test ~= 0
        
    %% loop through productivity states
    for j = 1:N
        
        % calculate  value at each level of current capital, for each level
        % of next period capital
        cons(:,:,j) = z(j).*kgrid.^alpha + (1-delta)*ones(n,1)*kgrid -...
         kgrid'*ones(1,n);
        util(:,:,j) = (cons(:,:,j).^(1-gamma)-1)/(1-gamma);

        % give large negative values to infeasible consumption choices
        % (i.e. enforce non-negativity constraint)
        utilj = util(:,:,j);
        
        Aj = cons(:,:,j);                % give negative value to 
        i = find( Aj <= 0);              % infeasible consumption choice
        utilj(i) = -10000;
        util(:,:,j) = utilj;
      
        % sum the value function at the current capital level
        futval = sum(v,2);
        futval = sum((ones(n,1)*P(j,:)).*v,2); % nkapx1 vector of
                                               % future value function

        vint(:,:,j) = util(:,:,j) + beta*futval*ones(1,n);
        [t1,t2] = max(vint(:,:,j));

        tv(j,:)     = t1; % RHS of Bellman
        tdecis(j,:) = t2; % index of capital level that gave tv
        
    end % end productivity state loop
    
    % calculate error and set up next iteration
    tv1     = tv';
    tdecis1 = tdecis';
    test= max(any(tdecis1-decis));
    err = norm(tv1-v);
    nprin=10; % print every 10 iterations
    if iter==nprin*round(iter/nprin)
        disp(['Iteration = ', num2str(iter),... 
            ', k Error = ', num2str(max(max(abs(tdecis1-decis)))),...
            ', V Error = ', num2str(err)])
    end
    iter  = iter + 1;
    decis = tdecis1;
    v     = tv1;

end
toc

%% calculate policy functions
kap   = kgrid(3:end); % save full grid
decis_ind = decis;
decis = (decis(3:end,:)-1)*inckap;
[xz,ykap] = meshgrid(z,kap);
condecis = xz.*ykap.^alpha + (1-delta).*ykap - decis;
v = v(3:end,:);
indvec = t2;

%% retrieve deterministic policy functions
params       = struct();
params.beta  = beta;
params.delta = delta;
params.z     = 1;
params.gamma = gamma;
params.alpha = alpha; % Cobb-Douglas exponent on capital

[kpol,cpol,k_ss,c_ss] = deterministic(params,kgrid);

%% define steady-state capital and consumption equations
%% (e.g. delta c = 0)
% delk0 = z*kap.^alpha - delta*kap;
% delc0 = z*kap.^alpha + (1-delta)*kap - knext;
% k_ss = knext;
% c_ss = z*k_ss^alpha - delta*k_ss;

%%% (a) display policy functions

%% plot consumption policy functions
figure(1)
    plot(kap,cpol,kap,condecis(:,1),kap,condecis(:,2))
    title('Policy function: c(k,z)')
    xlabel('k_t'); ylabel('c_t')
    legend(   'z = 1 (deterministic)',...
              'z = 0.8 (stochastic)',...
              'z = 1.2 (stochastic)',...
              'location','Southeast')
    saveas(gcf,'figure2ai.png')

%% plot capital policy functions
figure(2)
    plot(kap,kpol,kap,decis(:,1),kap,decis(:,2))
    title("Policy function: k'(k,z)")
    xlabel('k_t'); ylabel('k_{t+1}')
    legend(   'z = 1 (deterministic)',...
              'z = 0.8 (stochastic)',...
              'z = 1.2 (stochastic)',...
              'location','Southeast')
    saveas(gcf,'figure2aii.png')
      
%%% (b) simulate time series

%% determine number of burns and simulations
nBurn   = 100000;
nsim    = 100;

%% initialize 'observed data' vectors
obs_kap = zeros(nsim,1); obs_con = zeros(nsim,1); obs_w = zeros(nsim,1);
obs_r = zeros(nsim,1); obs_inv = zeros(nsim,1);

%% choose initial capital level and productivity state
k_ind = floor(length(kap)*rand(1,1)); % for laziness, an index of initial capital is chosen
j = 1;

%% simulate time series
for t = 1:(nBurn + nsim)
    
    %% determine next period's capital holdings
    k       = kap(k_ind);
    knext   = decis(k_ind,j);
    
    %% calculate current-period consumption
    c = condecis(k_ind,j);
    
    %% determine next period's productivity
    t_draw = rand(1,1); % random number between 0 and 1
    if t_draw <= P(j,1); j = 1; else; j = 2; end
    
    %% update k_ind for next period's simulation
    k_ind   = decis_ind(k_ind,j); %% knext index
    
    %% if past the burn phase, output time series
    if t > nBurn
        tt = t - nBurn;
        obs_kap(tt) = k;
        obs_con(tt) = c;
        obs_inv(tt) = z(j)*k^alpha-c;
        obs_w(tt)   = z(j)*(1-alpha)*k^alpha;
        obs_r(tt)   = z(j)*alpha*k^(alpha-1);
        fprintf('Simulation %d of %d\n',tt,nsim);
    else
        fprintf('Burning...\n')
    end
    
end


%% output results to .tex file
if exist('2b.tex', 'file')==2
  delete('2b.tex');
end
file2b = fopen('2b.tex','w');

tbl = ...
    ['\\begin{center}\n\\begin{tabular}{r c c}\n',...
    '& Deterministic & Stochastic \\\\ \\hline\n',...
    'Capital  & %4.3f & %4.3f  \\\\ \n',...
    'Investment & %4.3f & %4.3f  \\\\ \n',...
    'Consumption & %4.3f & %4.3f  \\\\ \n',...
    'Interest rate & %4.3f & %4.3f \\\\ \n', ...
    'Wage & %4.3f & %4.3f \\\\ \\hline\n',...
    '\\end{tabular}\n\\end{center}'];
fprintf(file2b,tbl,...
    round(k_ss,3),round(mean(obs_kap),3),...
    round(k_ss^alpha-c_ss,3),round(mean(obs_inv),3),...
    round(c_ss,3),round(mean(obs_con),3),...
    round(alpha*k_ss^(alpha-1),3),round(mean(obs_r),3),...
    round((1-alpha)*k_ss^alpha,3),round(mean(obs_w),3));
fclose(file2b);


%%% (c) calculate and output various moments from the simulated model

%% volatility of consumption, investment, and output
sd_c = std(obs_con);
sd_k = std(obs_kap);
sd_y = std(obs_kap.^alpha);

%% correlations
corr_c_y = corr(obs_con,obs_kap.^alpha);
corr_r_y = corr(obs_r,obs_kap.^alpha);

%% output results to .tex file
if exist('2c.tex', 'file')==2
  delete('2c.tex');
end
file2c = fopen('2c.tex','w');

tbl = ...
    ['\\begin{center}\n\\begin{tabular}{r c}\n',...
    '& \\\\ \\hline\n',...
    '$\\sigma_c$        & %4.3f \\\\ \n',...
    '$\\sigma_k$        & %4.3f \\\\ \n',...
    '$\\sigma_y$        & %4.3f \\\\ \n',...
    '$\\sigma_{c,y}$    & %4.3f \\\\ \n',...
    '$\\sigma_{r,y}$    & %4.3f \\\\ \\hline\n',...
    '\\end{tabular}\n\\end{center}'];
fprintf(file2c,tbl,...
    round(sd_c,3),round(sd_k,3),round(sd_y,3),...
    round(corr_c_y,3),round(corr_r_y,3));
fclose(file2c);


