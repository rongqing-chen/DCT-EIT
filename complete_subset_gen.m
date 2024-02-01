function [new_elem_centers, subset, coefficients_matrix] = complete_subset_gen(elem_centers, max_coeffs, basis_type)
%[new_elem_centers, subset, coefficients_matrix] = complete_subset_gen(
%   elem_centers, max_coeffs, basis_type) calculate the subset expansion 
%   and the corresponding new_elem_centers for the basis_type 
%   ('dct', 'poly'=chebyshev, 'cyl'). 
%   max_coeffs gives the maximum number of coefficient for every expansion.
%   coefficients_matrix is the matrix containing in every row the
%   combination of coefficients corresponding to the subset such that the
%   subset(:,i) corresponds to the coefficient_matrix(i)
%   
%   example: basis_type = {'dct', 'dct', 'poly'}
%   basis_type = {'cyl', ...} all the other entries are ignored.
 
[number_elements, number_dim] = size(elem_centers);

subset_makers = cell(number_dim,1);

new_elem_centers = zeros(size(elem_centers));

for idx = 1:length(basis_type)
    switch basis_type{idx}
        case 'dct'
            new_mins = pi/(2*number_elements);
            new_maxs = pi*(1-1/(2*length(elem_centers)));
            subset_makers{idx} = @make_DCT_subset;
            new_elem_centers(:,idx) = shift_elem_centers(elem_centers(:,idx), new_mins, new_maxs);
        case 'poly'
            new_mins = -1;
            new_maxs = +1;
            subset_makers{idx} = @make_chebyshev_subset;
            new_elem_centers(:,idx) = shift_elem_centers(elem_centers(:,idx), new_mins, new_maxs);
        case 'cyl'
            subset_makers = {@make_chebyshev_subset, @make_fourier_subset, @make_chebyshev_subset};
            new_elem_centers = cylindrical_elem_centers(elem_centers);
            new_elem_centers(:,1) = shift_elem_centers(new_elem_centers(:,1), -1, +1);
            new_elem_centers(:,3) = shift_elem_centers(new_elem_centers(:,3), -1, +1);
            break

        case 'cyl2' % this requires modifications
            subset_makers = {@make_chebyshev_subset, @make_fourier_subset, @make_chebyshev_subset};
            new_elem_centers = cylindrical_elem_centers(elem_centers);
            new_elem_centers(:,1) = shift_elem_centers(new_elem_centers(:,1), 0, +1);
            new_elem_centers(:,2) = shift_elem_centers(new_elem_centers(:,1), 0, pi);
            new_elem_centers(:,3) = shift_elem_centers(new_elem_centers(:,3), -1, +1);
            break
    end
end


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