function DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
%DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix)
% returns the nD rectangular DCT_subset calculated on elem_centers from the
% coefficients_matrix.

max_coeffs = max(coefficients_matrix);

DCT_subset = zeros(length(elem_centers), length(coefficients_matrix));

normalization_ = 2/(prod(max_coeffs));

for idk = 1:length(coefficients_matrix)
    a_jk = coefficients_matrix(idk,:);
    
    zero_freqs = sum(a_jk == 1);
    % each zero freq gets a 1/sqrt(2)
    normalization = normalization_/(sqrt(2).^zero_freqs); 

    DCT_subset(:,idk) = normalization*cos(a_jk(1)*elem_centers(:,1)).*cos(a_jk(2)*elem_centers(:,2));
end

end