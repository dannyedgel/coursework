%{
    This file is used to conduct all necessary tasks for problem set 3 of
    the second quarter of Econ 712

    Date created:  01 Dec 2020
    Last modified: 01 Dec 2020
    Author: Danny Edgel
%}

%%% clean workspace
clc; clear


%{
                    Question 3
   -------------------------------------------------
%}

%%% (b) solve for stationary distribution of l_t

%% define transition matrix
Q = [0.85, 0.15; 0.05, 0.95];

%% pick arbitrary starting point for l
l0 = [0, 1];

%% iterate through transitions until the probability of each labor state
%% changes by some very small amount (tolerance level)
tol = 1e-8;
i = 0; err = 1;
while (err > tol && i < 1000)
    i = i + 1;
    l = l0*Q;
    err = norm(l-l0);
    l0 = l;
    fprintf('Iteration %d, err = %d\n',i,err)
end

%% once finished, print stationary distribution of l
l

%%% (c) solve for the value function

%% define parameters
beta    = 0.95;
gamma   = 3;
r       = 0.03;
w       = 1.1;
N = 2;          % # of wage states
s=[.7,1.1];       % employment states

phi  = 0;  % borrowing limit,                               
                                                              
%% form capital grid
   
maxkap = 3;                      % maximum value of capital grid  
minkap = -phi;                   % borrowing constraint
inckap = 0.005;                   % size of capital grid increments
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

tic
while err > tol || test ~= 0 


    %  tabulate the utility function such that for zero or negative
    %  consumption utility remains a large negative number so that
    %  such values will never be chosen as utility maximizing      
    for j=1:N

      cons(:,:,j) = s(j)*w + (1+r)*ones(nkap,1)*kap -...
         kap'*ones(1,nkap);


      util(:,:,j) = (cons(:,:,j).^(1-gamma))/(1-gamma);
      utilj = util(:,:,j);

      Aj = cons(:,:,j);                % give negative value to 
      i = find(Aj <= 0);              % infeasible consumption choice
      utilj(i) = -10000;
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
   err = norm(tv1-v);
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
condecis = w*xs+(1+r)*ykap-decis;
kap=kap(3:end);
v = v(3:end,:);
indvec = t2;

% find abar
abar = min(kap(decis(:,2)<kap'));


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
    xlabel('a_t');ylabel('a_{t+1}'); xline(abar,'r-')
    text(abar-0.4,2.5,['$\overline{a}=$',num2str(abar)], ...
        'color','r','Interpreter', 'LaTeX')
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












