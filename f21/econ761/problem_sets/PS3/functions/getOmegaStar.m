% This function generates the OmegaStar matrix
%
% Inputs:
% x: matrix with firm IDs in the first column and brands owned by the firm
%       in the second column
% br_id: vector of brand IDs in the full data set
%
% Date created:  26 Oct 2021
% Last modified: 26 Oct 2021
% Author: Danny Edgel
%

function y = getOmegaStar(x, br_id)

brand_ind  = zeros(length(unique(x(:,2))));
uq_firm_id = unique(x(:, 1));
uq_brands  = unique(x(:, 2));

for i = 1:size(brand_ind, 1)
    for j = i:size(brand_ind, 1)
        for f = 1:length(uq_firm_id)
            rows = x(:, 1) == uq_firm_id(f);
            if1 = ismember(uq_brands(i), x(rows, 2));
            if2 = ismember(uq_brands(j), x(rows, 2));
            
            if if1 && if2; brand_ind(i, j) = 1; brand_ind(j, i) = 1; end
        end
    end
end

index = 1:size(brand_ind, 1);
y = zeros(length(br_id));

for i = 1:length(br_id)
    ibrand_ndx = index(uq_brands == br_id(i));
    for j = 1:length(br_id)
        jbrand_ndx = index(uq_brands == br_id(j));
        
        y(i, j) = brand_ind(ibrand_ndx, jbrand_ndx);
    end
end

end
