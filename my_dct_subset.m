% function (elem_centers)

% scale dimensions 
elem_centers = interp_mesh(fmdlStretch, 0); % center of elements
mins = min(elem_centers);
maxs = max(elem_centers);


M = 9;
N = 9;
norm_centers = (elem_centers - mins)./(maxs - mins);

orig_mins = [8, 38.72];
orig_maxs = [248, 217.28];
strecthed_centers = norm_centers.*(orig_maxs-orig_mins) + orig_mins;
new_centers(:,1) = pi*(2*strecthed_centers(:,1)+1)/(2*256);
new_centers(:,2) = pi*(2*strecthed_centers(:,2)+1)/(2*256);


values = zeros(length(new_centers), M*N);

idx = 0;
for jj = 0:M-1
    for kk = 0:N-1
        idx = idx+1;
        a_jk = [jj,kk];
        values(:,idx) = cos(a_jk(1)*new_centers(:,2)).*cos(a_jk(2)*new_centers(:,1));
    end
end

% values = sqrt(2/M)*sqrt(2/N)*values;
% 
% values(:,1) = values(:,1);

coeffs_2_check = [4,15];

imgRec.elem_data = values(:,coeffs_2_check);

close all
figure(1)
show_slices(imgRec)

imgRec.elem_data = spec_Mtx_col(:,coeffs_2_check);


figure(2)
show_slices(imgRec)
