function [norm_centers] = recenter_model(fmdl)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% scale dimensions 
fmdl_stretch = fmdl;
fmdl_stretch.nodes = fmdl.nodes;

% get position of fem elements
elem_centers = interp_mesh(fmdl_stretch, 0); % center of elements
mins = min(elem_centers);
maxs = max(elem_centers);

% normalize to 0-1 in both dimensions
norm_centers = (elem_centers - mins)./(maxs - mins);

end