function subset = complete_subset_gen(elem_centers, max_coeffs, basis_type)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[number_elements, number_dim] = size(elem_centers);

subset_makers = cell(number_dim,1);

new_mins = zeros(1,number_dim);
new_maxs = zeros(1,number_dim);

for idx = 1:length(basis_type)
    switch basis_type{idx}
        case 'dct'
            new_mins(idx) = pi/(2*number_elements);
            new_maxs(idx) = pi*(1-1/(2*length(elem_centers)));
            subset_makers{idx} = @make_DCT_subset;
        case 'poly'
            new_mins(idx) = -1;
            new_maxs(idx) = +1;
            subset_makers{idx} = @make_chebyshev_subset;
    end
end

[new_elem_centers] = shift_elem_centers(elem_centers, new_mins, new_maxs);

% coefficients ordered in row, by col, by depth
if number_dim == 2
    [MM, NN] = ndgrid(1:max_coeffs(1), 1:max_coeffs(2));
    coefficients_matrix = [MM(:), NN(:)];
elseif number_dim == 3
    % coefficients ordered in row, by col, by depth
    [MM, NN, OO] = ndgrid(1:max_coeffs(1), 1:max_coeffs(2), 1:max_coeffs(3));
    coefficients_matrix = [MM(:), NN(:), OO(:)];
end

subset = make_subset_3D(new_elem_centers, coefficients_matrix, subset_makers);

end