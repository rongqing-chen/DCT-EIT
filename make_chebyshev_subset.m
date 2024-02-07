function poly_subset = make_chebyshev_subset(elem_centers, coefficients_matrix)
%poly_subset = make_chebyshev_subset(elem_centers, coefficients_matrix)
% make a nD rectangular chebyshev subset from the elem_centers (centers of
% fem elements)
% and the a coefficient matrix. 

max_coeffs = max(coefficients_matrix, [],1);

[combinations_of_coeffs, number_of_dimensions] = size(coefficients_matrix);
single_dim_subsets = cell(number_of_dimensions,1);
poly_subset = ones(length(elem_centers), combinations_of_coeffs);

for idm = 1:number_of_dimensions
    single_dim_subsets{idm} = chebyshev_polynomials(elem_centers(:,idm), max_coeffs(idm));
end

for idk = 1:combinations_of_coeffs
    for idm = 1:number_of_dimensions
        poly_subset(:,idk) = poly_subset(:,idk).*single_dim_subsets{idm}(:,coefficients_matrix(idk,idm));
    end
end

end