%{
    This file performs all analayses requested in part 1 of problem set 3
    for Econ 761

    Date created:  23 Oct 2021
    Last modified: 24 Oct 2021
    Author: Danny Edgel
%}
clear; clc

%% Load the data and functions
load('data/ps2.mat'); load('data/iv.mat')

addpath('downloaded/code')
addpath('functions')


%% Prepare data for MNL estimation

% normalize by first brand
norm_id = floor(id(1)/100000);
new_id  = id((floor(id/100000)~=norm_id));
cq      = (id - (floor(id/100000)*100000));
uq_cq   = unique(cq);
new_cq  = (new_id - (floor(new_id/100000)*100000));
X = full([x1(floor(id/100000) ~= norm_id, 1), ...
    x1(floor(id/100000) ~= norm_id, 3:end)]);
Z = iv(floor(id/100000) ~= norm_id, 2:end);
S = zeros(size(X, 1), 1);

for i = 1:length(uq_cq)
    new_i   = new_cq == uq_cq(i);
    old_i   = cq == uq_cq(i) & floor(id/100000) ~= norm_id;
    norm_i  = cq == uq_cq(i) & floor(id/100000) == norm_id;
    X(new_i, 1) = X(new_i, 1) - x1(norm_i, 1);
    S(new_i, 1) = log(s_jt(old_i)) - log(s_jt(norm_i));
end

fbr_id    = floor(id/100000);
firm_id   = floor(fbr_id/1000);
uq_firms  = unique(firm_id);
yqcity    = id - fbr_id * 100000;
uq_yqcity = unique(yqcity);

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
uq_brands  = unique(fbr_id - firm_id*1000);
uq_fbr     = unique(fbr_id);
uq_firm_id = unique(firm_id);
br_id      = fbr_id - firm_id*1000;
dat        = [floor(uq_fbr/1000), uq_fbr - floor(uq_fbr/1000)*1000];

index = 1:length(uq_brands);
OmegaStar = getOmegaStar(dat, br_id);

% generate matrices with repeated share vectors for easy computing in-loop
s_rows = repmat(s_jt, 1, length(s_jt)); s_cols = s_rows';
p = full(x1(:, 1));

for i = 1:length(alphas)
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




