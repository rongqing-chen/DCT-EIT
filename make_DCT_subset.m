function DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
%DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
% returns the nD rectangular DCT_subset calculated on elem_centers from the
% coefficients_matrix. elem_centers should between 0 and pi.

max_coeffs = max(coefficients_matrix);

[number_of_coeffs, number_of_dimensions] = size(coefficients_matrix);
DCT_subset = zeros(length(elem_centers), number_of_coeffs);

normalization_ = sqrt(2).^number_of_dimensions/prod(max_coeffs);

for idk = 1:number_of_coeffs
    a_jk = coefficients_matrix(idk,:);
    
    zero_freqs = sum(a_jk == 1);
    % each zero freq gets a 1/sqrt(2)
    normalization = normalization_/(sqrt(2).^zero_freqs); 
    DCT_subset(:,idk) = normalization*prod(cos(elem_centers.*(a_jk-1)),2);
end

end