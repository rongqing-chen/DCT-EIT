function ordered_coeffs = order_coeffs_tensor_product(coeffs_dir_1, coeffs_dir_2)
%ordered_coeffs = order_coeffs_tensor_product(coeffs_dir_1, coeffs_dir_2)
% return the ordered_coeffs for a 2D rectangular tensor product with 
% coeffs_dir_1 and coeffs_dir_2. The ordered_coeffs are ordered according
% to the complessive order.

[MM, NN] = ndgrid(coeffs_dir_1, coeffs_dir_2);
MM = zigzag_matrix(MM);
NN = zigzag_matrix(NN);
ordered_coeffs = [MM(:), NN(:)];

end