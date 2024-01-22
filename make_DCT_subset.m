function DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
%DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
% returns the 2D rectangular DCT_subset calculated on elem_centers from the
% coefficients_matrix.

x_max_coeff = max(coefficients_matrix(:,1));
y_max_coeff = max(coefficients_matrix(:,2));
DCT_subset = zeros(length(elem_centers), length(coefficients_matrix));

for idk = 1:length(coefficients_matrix)
    a_jk = coefficients_matrix(idk,:);
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