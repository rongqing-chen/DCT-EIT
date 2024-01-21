function [DCT_subset, ordered_coefficients] = make_DCT_subset(elem_centers, x_max_coeff, y_max_coeff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% order coefficients in zig zag way
ordered_coefficients = order_coeffs_tensor_product(0:y_max_coeff-1, 0:y_max_coeff-1);

DCT_subset = zeros(length(elem_centers), x_max_coeff*y_max_coeff);

for idk = 1:length(ordered_coefficients)
    a_jk = ordered_coefficients(idk,:);
    normalization = 2/(x_max_coeff*y_max_coeff);
    if a_jk(1) == 0
        normalization = normalization/sqrt(2);
    end
    if a_jk(2) == 0
        normalization = normalization/sqrt(2);
    end
    DCT_subset(:,idk) = normalization*cos(a_jk(1)*elem_centers(:,1)).*cos(a_jk(2)*elem_centers(:,2));
end

end