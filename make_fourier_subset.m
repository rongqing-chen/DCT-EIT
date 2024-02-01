function subset = make_fourier_subset(elem_centers, coefficients_matrix)
%subset = make_fourier_subset(elem_centers, coefficients_matrix)
% make a 1D fourier subset from the elem_centers (centers of
% fem elements)
% and the a coefficient matrix. 

max_coeffs = max(coefficients_matrix, [],1);
if mod(max_coeffs,2) ~= 1
    error('The max coefficient (%d) was not odd', max_coeffs)
end

[combinations_of_coeffs, number_of_dimensions] = size(coefficients_matrix);

subset = zeros(length(elem_centers), combinations_of_coeffs);

cos_basis = cos(elem_centers.*(1:floor(max_coeffs/2)));
sin_basis = sin(elem_centers.*(1:floor(max_coeffs/2)));

total_basis = zeros(length(elem_centers), max_coeffs);
total_basis(:,1) = 1;
total_basis(:,2:2:end) = cos_basis;
total_basis(:,3:2:end) = sin_basis;

for idk = 1:combinations_of_coeffs
    subset(:,idk) = total_basis(:,coefficients_matrix(idk,1));
end

end