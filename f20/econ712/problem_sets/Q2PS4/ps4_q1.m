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
r       = 0.03;
w       = 1.1;
N = 2;          % # of wage states
s=[.7,1.1];       % employment states

phi  = 0;  % borrowing limit,    

%% define firm parameters
alpha = 0.36;
delta = 0.08;
                                                              
%% form capital grid
   
maxkap = 7;                      % maximum value of capital grid  
minkap = -phi;                   % borrowing constraint
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

%  initialize some variables
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
while err > tol || test ~= 0 


    %  tabulate the utility function such that for zero or negative
    %  consumption utility remains a large negative number so that
    %  such values will never be chosen as utility maximizing      
    for j=1:N
        
        % calculate wages and interest rates for each capital decision
        w = (1-alpha)*(kap/s(j)).^alpha;
        r = kap.^(alpha-1)*(s(j)^(1-alpha));
        
        cons(:,:,j) = s(j).*w + (1+r).*(ones(nkap,1)*kap) -...
            kap'*ones(1,nkap);
        util(:,:,j) = (cons(:,:,j).^(1-gamma))/(1-gamma);
        utilj = util(:,:,j);
        
        Aj = cons(:,:,j);
        utilj(Aj <= 0) = -10000;
        util(:,:,j) = utilj;
        
        futval = sum((ones(nkap,1)*Q(j,:)).*v,2); % nkap by 1 vector of
                                                  % future value function
        vint(:,:,j) = util(:,:,j) + beta*futval*ones(1,nkap);
        [t1,t2] = max(vint(:,:,j));
        tv(j,:) =t1;
        tdecis(j,:)=t2;
    end % end labor states loop
    
    tv1     = tv';
    tdecis1 = tdecis';
    test=max(any(tdecis1-decis));
    diff = tv1-v;
    err = norm(diff(~isnan(diff)));
    nprin=10;
    if iter==nprin*round(iter/nprin)
        disp(['Iteration=', num2str(iter),... 
            ', k Error=', num2str(max(max(abs(tdecis1-decis)))),...
            ', V Error=', num2str(err)])
    end
    iter=iter+1;
    decis=tdecis1;
    v=tv1;
end % end outermost while loop
toc

% focus on decisions off lowest boundary...
kgrid = kap; % save full grid
decis = -phi + (decis(3:end,:)-1)*inckap;
[xs,ykap] = meshgrid(s,kap(3:end));

% calculate wages and interest rates for each capital and labor level
w = (1-alpha)*(ykap./xs).^alpha;
r = ykap.^(alpha-1).*(xs.^(1-alpha));

condecis = w.*xs+(1+r).*ykap-decis;
kap=kap(3:end);
v = v(3:end,:);
indvec = t2;

%% Calculate equilibrium--simulate market


%% output results to .tex file
if exist('1a.tex', 'file')==2
  delete('1a.tex');
end
file1b = fopen('1a.tex','w');

tbl = ...
    ['\\begin{center}\n\\begin{tabular}{r c c}\n',...
    '&  \\\\ \\hline\n',...
    'Capital  & %4.3f \\\\ \n',...
    'Interest rate & %4.3f \\\\ \n', ...
    'Wage & %4.3f \\\\ \\hline\n',...
    '\\end{tabular}\n\\end{center}'];
fprintf(file1b,tbl,...
    round(k_ss,3),round(mean(obs_kap),3),...
    round(k_ss^alpha-c_ss,3),round(mean(obs_inv),3),...
    round(c_ss,3),round(mean(obs_con),3),...
    round(alpha*k_ss^(alpha-1),3),round(mean(obs_r),3),...
    round((1-alpha)*k_ss^alpha,3),round(mean(obs_w),3));
fclose(file1b);

% (i) plot value funcitons for high and low labor states
figure(1)
    plot(kap,v(:,1),kap,v(:,2))
    title('Value functions')
    xlabel('a_t'); ylabel('V(a,l)')
    legend('l_i=0.7','l_i=1.1','location','Southeast')
    saveas(gcf,'figure1.png')
    
% (ii) Plot asset decision rules for each labor state
figure(2)
    plot(kap,decis(:,1),kap,decis(:,2),kap,kap,'k--')
    title("Policy functions: a'(a,l_i)")
    xlabel('a_t');ylabel('a_{t+1}');
    legend('l_i=0.7','l_i=1.1','a_t=a_{t+1}','location','Southeast')
    saveas(gcf,'figure2.png')
   
%% (d) Find stationary distribution of asset holdings and plot

% initialize arbitrary asset distribution
adist0       = zeros(nkap,N);
adist0(1)    = 1;
err = 1; % initialize tolerance
i   = 0;
while (err > tol && i < 1000)
    adist = 0*adist0;
    
    % loop through labor states and starting asset amounts
    for n = 1:N
       for a = 1:nkap
           if adist0(a,n)>1e-14 % ignore empty asset holdings
               targ = tdecis1(a,n);
               adist(tdecis1(a,n),1) = adist(targ,1)+adist0(a,n)*Q(n,1);
               adist(tdecis1(a,n),2) = adist(targ,2)+adist0(a,n)*Q(n,2);
           end
       end
    end
    i = i+1;
    err = norm(adist0(:)-adist(:));
    adist0 = adist;
    fprintf('Iteration %d, err = %d\n',i,err)
end
amarg = sum(adist,2);
ave = sum(amarg.*kgrid');

% plot marginal asset distribution
figure(3)
    bar(kgrid,amarg)
    title('Marginal distribution of a_t')
    xlim([-0.1 max(kgrid)]); xline(ave,'r-')
    text(ave+0.1,0.06,['Mean=',num2str(ave)], 'color','r')
    xlabel('a_t'); ylabel('pmf')
    saveas(gcf,'figure3.png')

