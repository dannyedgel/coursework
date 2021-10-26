% This function calculates multi-product Bertrand-Nash markups
%
% Inputs:
% x: vector of brand IDs, sorted and indexed by quality, owned by the same
%       firm
% a: the alpha value from an MNL regression
% s: vector of market shares for the firm for each product
% i: start index of x
% mu: also the output; vector of multi-product Bertrand-Nash markups
%
% Date created:  26 Oct 2021
% Last modified: 26 Oct 2021
% Author: Danny Edgel
%

function mu = BertNashMu(x, a, s, i, mu)
    if x(i) ~= x(i + 1) - 1
        mu(i) = ((1-s(i))*(-a))^(-1);
        mu = BertNashMu(x, i + 1, mu);
    else
        j = i;
        k = 0;
        while k == 0 && j <= length(x)
           if  x(i) ~= x(i + 1) - 1; j = 
               
        end
    end

end