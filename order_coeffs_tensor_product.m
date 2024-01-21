function ordered_coeffs = order_coeffs_tensor_product(coeffs_dir_1, coeffs_dir_2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[MM, NN] = ndgrid(coeffs_dir_1, coeffs_dir_2);

ordered_coeffs = [zigzag_matrix(MM)', zigzag_matrix(NN)'];

% end