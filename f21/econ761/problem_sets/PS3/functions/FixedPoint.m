% This function solves for the fixed point of the post-merger price vector
%
% Inputs:
% Own: Post-merger ownership matrix
% H: Original matrix of cross-partials
% mc: vector of estimated marginal costs
% p0: initial price vector guess
% s0: pre-merger shares
% a: alpha hat
%
% Outputs:
% p: fixed point of price vector
% s: shares associated with fixed point of price vector
% Omega: Omega matrix associated with fixed point of price vector
%
% Date created:  27 Oct 2021
% Last modified: 27 Oct 2021
% Author: Danny Edgel
%

function [p, s, Omega] = FixedPoint(own, H, mc, p0, s0, a)

% initialize oshares for first guess
Omega = own .* H; s = s0;
err = 100;

% iteratively solve for shares and re-calculate prices until convergence
i = 1;
while err > 10-8
    
    % update Omega according to last loop's shares
    s_rows = repmat(s, 1, length(s)); s_cols = s_rows';
    
    H = (-a*(s_rows.*s_cols)).*ones(size(own));
    H = H - diag(diag(H)) + diag(s.*(1-s)*a);
    Omega = own .* H;
    
    % update p and s
    p = mc + Omega\s;
    s = Omega*(p0 - mc);
    
    err = norm(p - p0); p0 = p;
    
    %disp(['Iteration ', string(i), ', err = ', string(err)])
    i = i + 1;
end

end