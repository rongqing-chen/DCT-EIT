function unstruct_mask = make_unstructured_mask(fmdl, mask_matrix)
%unstruct_maks = make_unstructured_mask(fmdl, mask_matrix)
% map a rectangular mask_matrix into fmdl and return unstruct_mask the 
% masked points on the elements. 

pts = interp_mesh(fmdl, 0);

[x_dim, y_dim] = size(mask_matrix);
x = 0:x_dim-1;
y = (0:y_dim-1)';

mask_interp = griddedInterpolant({x,y}, mask_matrix, 'nearest');

% dimensions need inversion
unstruct_mask = mask_interp([pts(:,2), pts(:,1)]);


end