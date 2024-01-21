function [fmdlStretch, new_centers] = scale_model_dimension(fmdl, magic_values)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% magic values
orig_mins = [8, 38.72]; % magic_values(1,:);
orig_maxs = [248, 217.28]; % magic_values(2,:);

% scale dimensions 
fmdlStretch = fmdl;
fmdlStretch.nodes = fmdl.nodes * 120 + 256/2; % fmdl.nodes * magic_values(3,1) + magic_values(3,2);

% get position of fem elements
elem_centers = interp_mesh(fmdlStretch, 0); % center of elements
mins = min(elem_centers);
maxs = max(elem_centers);

% normalize to 0-1 in both dimensions
norm_centers = (elem_centers - mins)./(maxs - mins);

% remap to new positions
strecthed_centers = norm_centers.*(orig_maxs-orig_mins) + orig_mins;
new_centers(:,1) = pi*(2*strecthed_centers(:,1)+1)/(2*256);
new_centers(:,2) = pi*(2*strecthed_centers(:,2)+1)/(2*256);


end