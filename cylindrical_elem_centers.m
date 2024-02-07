function cyl_elem_centers = cylindrical_elem_centers(elem_centers)
%cyl_elem_centers = cylindrical_elem_centers(elem_centers) gives the
%   cyl_elem_centers in cylindrical coordinates (r, angle, height). 
%   It works for 2 and 3D

n_dimensions = size(elem_centers, 2);

[new_elem_centers] = shift_elem_centers(elem_centers(:,1:2), [-1,-1], [1,1]);

[angle_coords, radius_coords] = cart2pol(new_elem_centers(:,1), new_elem_centers(:,2));

if n_dimensions == 2
    cyl_elem_centers = [radius_coords, angle_coords];
elseif n_dimensions == 3
    cyl_elem_centers = [radius_coords, angle_coords, elem_centers(:,3)];
end

end