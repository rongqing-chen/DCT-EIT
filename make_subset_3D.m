function [subset] = make_subset_3D(elem_centers, coefficients_matrix, subset_makers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

single_dim_subsets = zeros(length(elem_centers), size(coefficients_matrix,1), 3);
for idm = 1:3
    single_dim_subsets(:,:,idm) = subset_makers{idm}(elem_centers(:,idm),coefficients_matrix(:,idm));
end

subset = prod(single_dim_subsets,3);

end