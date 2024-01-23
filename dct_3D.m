


elem_centers = recenter_model(fmdl)*pi;

noise_ampl = 0.01;
vi_noise = vi;
vi_noise.meas = vi.meas + noise_ampl*std(vi.meas).*randn(size(vi.meas)); 
dva = calc_difference_data( vh, vi_noise, fmdl);

M = 12;
N = 12;
O = 12;

% coefficients ordered in row, by col, by depth
[MM, NN, OO] = ndgrid(0:M-1, 0:N-1, 0:O-1);
coefficients_matrix = [MM(:), NN(:), OO(:)];

DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix);

J_DCT = J* DCT_subset;

R = eye(size(J_DCT,2));
%%
lambda = 1e-4; %without noise you can use 1e-6 or less
dctCoeff = (J_DCT'*J_DCT + lambda.^2*R)\(J_DCT'*dva);

%% inverse DCT
recCond = DCT_subset*dctCoeff;

%% assigning element values

imgr.elem_data = recCond;


figure(4)
subplot(1,2,1)
show_fem(imgr);


subplot(1,2,2)
show_3d_slices(imgr, [1,1.9], [0.5],[0.5]);
view(-14,13); axis tight; axis equal;


