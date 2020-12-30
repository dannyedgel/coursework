
%%% parameters
syms a b f


%%% social planner's allocation
t_sp = 1/2;
ell_sp = (2/3)*(1-a);
l_sp = 1-ell_sp;
c_sp = (1-t_sp)*ell_sp;
g_sp = t_sp*ell_sp;

u_sp = log(l_sp) + log(a + c_sp) + log(a + g_sp);

%%% RE allocation
t_r = (1/2)*(4-3*a-sqrt(12-20*a+9*(a^2)));
ell_r = (1-t_r-a)/(2*(1-t_r));
l_r = 1-ell_r;
c_r = (1/2)*(1-t_r-a);
g_r = t_r*ell_r;

u_r = log(l_r) + log(a + c_r) + log(a + g_r);

%% NE allocation
t_n = (1/2);
l_n = (1/2)*(1+2*a);
ell_n = 1- l_n;
c_n = (1-t_n)*ell_n;
g_n = t_n*ell_n;

u_n = log(l_n) + log(a + c_n) + log(a + g_n);

%%% RE deviation
ell_r = (1-t_r-a)/(2*(1-t_r));
l_r = 1-ell_r;
c_r = (1/2)*(1-t_sp-a);
g_r = t_sp*ell_r;

u_dev = log(l_r) + log(a + c_r) + log(a + g_r);


%%% 

%%% Test inequality
eqn = (b/(1-b))*(u_r-u_n)-(u_dev-u_r)==f;

[S,param,cond] = solve(eqn,f>0,'ReturnConditions',true);

alpha = linspace(0,0.5,10000);
alpha = alpha(2:(end-1));

figure(1)
plot1 = plot(alpha,subs(S,'a',alpha,'b',.99));