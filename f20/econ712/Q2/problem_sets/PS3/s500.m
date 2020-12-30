
%% Example of discrete state dynamic programming
%% adapted from example in Ljungvist and Sargent 
%%    and code from Sargent's web page


mu     = 3; % risk aversion \\
beta   = 0.9;   % subjective discount factor \\
wage= .2; r=0.05; % wages and interest rates \\
N = 2; prob=[.8,.2;.05,.95]; % \# of wage states and transitions\\
s=[.2,1]; % employment states

      phi  = (wage*s(1)/r);          
% -phi is borrowing limit, 
%--------------------------------------------------------------------                               
                                                              
%   form capital grid
   
maxkap = 3;                     % maximum value of capital grid  
minkap = -phi;                   % borrowing constraint
inckap = 0.02;                    % size of capital grid increments
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

   %  initialize some variables
   %
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
   while test ~= 0;
      
  
   %  tabulate the utility function such that for zero or negative
   %  consumption utility remains a large negative number so that
   %  such values will never be chosen as utility maximizing      
   for j=1:N
      
      
      cons(:,:,j) = s(j)*wage + (1+r)*ones(nkap,1)*kap -...
         kap'*ones(1,nkap);
            
            
      util(:,:,j) = (cons(:,:,j).^(1-mu)-1)/(1-mu);
      utilj = util(:,:,j);
      
      Aj = cons(:,:,j);                % give negative value to 
      i = find( Aj <= 0);              % infeasible consumption choice
      utilj(i) = -10000;
      util(:,:,j) = utilj;

      futval = sum((ones(nkap,1)*prob(j,:)).*v,2); % nkap by 1 vector of
                                                   % future value function
                                                   
                                             
      vint(:,:,j) = util(:,:,j) + beta*futval*ones(1,nkap);
      [t1,t2] = max(vint(:,:,j));
      
      tv(j,:) =t1;
      tdecis(j,:)=t2;
   end
       tv1     = tv';
       tdecis1 = tdecis';
       test=max(any(tdecis1-decis));
       nprin=10;
       if iter==nprin*round(iter/nprin)
           disp(['Iteration=', num2str(iter),... 
               ', k Error=', num2str(max(max(abs(tdecis1-decis)))),...
               ', V Error=', num2str(max(max(abs(tv1-v))))])
       end
       iter=iter+1;
       decis=tdecis1;
       v=tv1;
   end;
   toc
   % focus on decisions off lowest boundary...
   
   decis = -phi + (decis(3:end,:)-1)*inckap;
   [xs,ykap] = meshgrid(s,kap(3:end));
   condecis = wage*xs+(1+r)*ykap-decis;
   kap=kap(3:end);
  
figure(1)
      plot(kap,decis(:,1),kap,decis(:,2),kap,kap,'k--')
      title('Policy function: Capital')
      xlabel('asset of current period')
      ylabel('asset of next period')
      legend('Low s','High s','a_t=a_{t+1}','location','Northwest')
      print -depsc pol_cap.eps
      
figure(2)
      plot(kap,condecis(:,1),kap,condecis(:,2))
      title('Policy function: Consumption')
      xlabel('current asset')
      ylabel('consumption')
      legend('Low s','High s','location','Northwest')
      print -depsc pol_cons.eps
   

