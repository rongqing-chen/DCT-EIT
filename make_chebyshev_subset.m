function poly_subset = make_chebyshev_subset(elem_centers, coefficients_matrix)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

x_max_coeff = max(coefficients_matrix(:,1));
y_max_coeff = max(coefficients_matrix(:,2));
poly_subset = zeros(length(elem_centers), length(coefficients_matrix));

bx =  chebyshev_polynomials(elem_centers(:,1), x_max_coeff+1);
by =  chebyshev_polynomials(elem_centers(:,2), y_max_coeff+1);


for idk = 1:length(coefficients_matrix)
    a_jk = coefficients_matrix(idk,:);
    poly_subset(:,idk) = bx(:,a_jk(1)+1).*by(:,a_jk(2)+1);
end

end