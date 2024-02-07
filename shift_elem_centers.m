function [new_elem_centers] = shift_elem_centers(elem_centers, new_mins, new_maxs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

mins = min(elem_centers);
maxs = max(elem_centers);

% normalize to 0-1 in both dimensions
norm_centers = (elem_centers - mins)./(maxs - mins);

new_elem_centers = norm_centers.*(new_maxs - new_mins) + new_mins;

end