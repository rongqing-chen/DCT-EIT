function unstruct_maks = make_unstructured_mask(fmdl, mask_matrix)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

pts = interp_mesh(fmdl, 0);

[x_dim, y_dim] = size(mask_matrix);
x = 0:x_dim-1;
y = (0:y_dim-1)';

mask_interp = griddedInterpolant({x,y}, mask_matrix, 'nearest');

% dimensions need inversion
unstruct_maks = mask_interp([pts(:,2), pts(:,1)]);


end