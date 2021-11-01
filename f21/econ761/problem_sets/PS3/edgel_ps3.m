%{
    This file performs all analayses requested in part 1 of problem set 3
    for Econ 761

    Date created:  23 Oct 2021
    Last modified: 31 Oct 2021
    Author: Danny Edgel
%}
clear; clc

%% Load the data and functions
load('data/ps2.mat'); load('data/iv.mat')

addpath('downloaded/code')
addpath('functions')

global invA ns x1 x2 s_jt IV vfull dfull theta1 theti thetj cdid cdindex


%% Prepare data for MNL estimation
N       = length(id);
norm_id = floor(id(1)/100000);
cq      = (id - (floor(id/100000)*100000));
uq_cq   = unique(cq);

fbr_id    = floor(id/100000);
firm_id   = floor(fbr_id/1000);
uq_firms  = unique(firm_id);
yqcity    = id - fbr_id * 100000;
uq_yqcity = unique(yqcity);

uq_brands  = unique(fbr_id - firm_id*1000);
uq_fbr     = unique(fbr_id);
uq_firm_id = unique(firm_id);
br_id      = fbr_id - firm_id*1000;

nbr = length(uq_fbr); ncq = length(uq_cq);

cid = floor(cq/1000); cdid = ones(length(cid), 1);
for i = 2:N; if cid(i) ~= cid(i-1); cdid(i) = cdid(i-1)+1; 
    else; cdid(i) = cdid(i-1); end; end

% normalize by "outside option"--the market share not captured by in-sample
% firms (note: this code taken from Nevo)
temp = cumsum(s_jt);
cdindex = (nbr:nbr:nbr*ncq)';
sum1 = temp(cdindex,:);
sum1(2:size(sum1,1),:) = diff(sum1);
outshr = 1.0 - sum1(ones(length(floor(cq/1000)), 1),:);


X = full(x1);
Z = [iv, X(:, 2:end)];
S = log(s_jt) - log(outshr);

%% Estimate specifications

% OLS
est.ols = ols(S, X(:, 1), 1);

% OLS with brand FE
est.fe  = ols(S, X, 1);

% IV without brand FE
est.iv = regressIV(S, X(:, 1), [], Z, 1);

% IV with brand FE
est.iv_fe = regressIV(S, X(:, 1), X(:, 3:end), Z, 1);



% Print results
printTable1(est, S);


%% Compute markups
markups = struct();
names = {'ols', 'fe', 'iv', 'iv_fe'};
alphas = [est.ols.b(1), est.fe.b(1), est.iv.b(1), est.iv_fe.b(1)];

%OmegaStar = zeros(length();

% generate JxJ matrix where each entry =1 if j and k are owned by the same
% firm and =0 otherwith
dat        = [floor(uq_fbr/1000), uq_fbr - floor(uq_fbr/1000)*1000];

index = 1:length(uq_brands);
OmegaStar = getOmegaStar(dat, br_id);

% generate matrices with repeated share vectors for easy computing in-loop
s_rows = repmat(s_jt, 1, length(s_jt)); s_cols = s_rows';
p = full(x1(:, 1));

for i = 1:length(alphas)
   markups.(names{i})       = est.(names{i});
   markups.(names{i}).alpha = alphas(i);
   
   % generate a JxJ matrix of cross-price elasticities for current alpha
   H = (-alphas(i)*(s_rows.*s_cols)).*ones(size(OmegaStar));
   H = H - diag(diag(H)) + diag(s_jt.*(1-s_jt)*(alphas(i)));
   
   Omega = OmegaStar.*H;
   
   % calculate markups according to Lecture 9, slide 23:
   markups.(names{i}).mu = Omega\s_jt;
   
   % back out marginal costs and margins
   markups.(names{i}).mc        = p - markups.(names{i}).mu;
   markups.(names{i}).margins   = p ./ markups.(names{i}).mc - 1;
end


% Print results
printTable2(markups, names);


%% Simulate post-merger equilibrium

% prepare results struct
mergeNames = {'PostNabisco', 'GMQuaker'}; mergeType = {[3, 6], [2, 4]};
mergers.None = markups; mergers.None.OmegaStar = OmegaStar;
mergers.None.Omega = Omega;

% generate new ownership structures for each merger scenario
for m = 1:length(mergeNames)
    mergers.(mergeNames{m}) = markups;
    
    % make merger adjustment in ownership matrix
    x = dat;
    x(dat(:, 1) == mergeType{m}(2), 1) = mergeType{m}(1);
    OmegaStar = getOmegaStar(x, br_id);

    % calculate mu and price for each alpha estimate
    for i = 1:length(alphas)
       H = (-alphas(i)*(s_rows.*s_cols)).*ones(size(OmegaStar));
       H = H - diag(diag(H)) + diag(s_jt.*(1-s_jt)*(alphas(i)));
       
       [p, s, Omega] = FixedPoint(OmegaStar, H, ...
           mergers.(mergeNames{m}).(names{i}).mc, p, s_jt, alphas(i));


       % calculate markups, margins, etc:
       mergers.(mergeNames{m}).(names{i}).mu = ...
           p - mergers.(mergeNames{m}).(names{i}).mc;
       mergers.(mergeNames{m}).(names{i}).p = p;
       mergers.(mergeNames{m}).(names{i}).margins = ...
           mergers.(mergeNames{m}).(names{i}).p ./ ...
           mergers.(mergeNames{m}).(names{i}).mc - 1;
       mergers.(mergeNames{m}).(names{i}).Omega = Omega;
       mergers.(mergeNames{m}).(names{i}).s = s;
    end
    
    % save simulation results in merger struct
    mergers.(mergeNames{m}).OmegaStar = OmegaStar;
end

% Print results
printTable3(mergers, names, mergeNames);

%% Estimate mixed logit model

% Use Nevo's BLP code, with some edits to adapt it to my syntax and for
% efficiency

ns = 20;       % number of simulated "individuals" per market % 
IV = Z;

% input starting values for non-linear variable coefficients
theta2w=    [0.3772    3.0888         0    1.1859         0;
             1.8480   16.5980    -.6590         0   11.6245;
            -0.0035   -0.1925         0    0.0296         0;
             0.0810    1.4684         0   -1.5143         0];
    
% create a vector of the non-zero elements in the above matrix, and the %
% corresponding row and column indices. this facilitates passing values % 
% to the functions below. %
[theti, thetj, theta2]=find(theta2w);


horz=['    mean       sigma      income   income^2    age    child'];
vert=['constant  ';
      'price     ';
          'sugar     ';
          'mushy     '];

% create weight matrix
A       = Z'*Z;
invA    = inv(A);

% Logit results and save the mean utility as initial values for the search below

mid = x1'*Z*(A\Z');
t = ((mid*x1)\mid)*S;        % IV of log shares on X1 using IV as instruments
mvalold = x1*t;              % Fitted log shares
oldt2 = zeros(size(theta2)); % Zero out old theta2
mvalold = exp(mvalold);      % Compute shares

save mvalold mvalold oldt2
clear mid y outshr t oldt2 mvalold temp sum1

vfull = v(cdid,:);
dfull = demogr(cdid,:);

options = optimset('GradObj','on','TolFun',0.1,'TolX',0.01);

% the following line computes the estimates using a Quasi-Newton method % 
% with an *analytic* gradient %
[theta2,fval,exitflag,output] = fminunc('gmmobjg',theta2,options);

% computing the s.e.
vcov = var_cov(theta2);
se = sqrt(diag(vcov));

theta2w = full(sparse(theti,thetj,theta2));
t = size(se,1) - size(theta2,1);
se2w = full(sparse(theti,thetj,se(t+1:size(se,1))));

% print results
preMerger.vcov = vcov; preMerger.theta1 = theta1; preMerger.x2 = x2;
preMerger.se = se; preMerger.theta2w = theta2w; preMerger.fval = fval;
preMerger.se2w = se2w;

printTable4(preMerger);




