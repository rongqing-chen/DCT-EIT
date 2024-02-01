clear
run reconstruction_3D

%% create second model (avoid inverse crime)
fmdl2 = ng_mk_cyl_models([3,1,0],[15,1,1.5,2],[0.1,0,0.05]); % 5293 nodes
imdl2 = imdl;
imdl2.fwd_model = fmdl2;
imdl2.fwd_model.stimulation = imdl.fwd_model.stimulation;
img_inv = mk_image(imdl2,1);

% calc the Jacobian
J = calc_jacobian(img_inv);

%%
elem_centers = recenter_model(fmdl2)*pi;

noise_ampl = 0.01;
vi_noise = vi;
vi_noise.meas = vi.meas + noise_ampl*std(vi.meas).*randn(size(vi.meas)); 
dva = calc_difference_data( vh, vi_noise, fmdl);

M = 12;
N = 12;
O = 12;

% coefficients ordered in row, by col, by depth
[MM, NN, OO] = ndgrid(1:M, 1:N, 1:O);
coefficients_matrix = [MM(:), NN(:), OO(:)];

DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix);

J_DCT = J* DCT_subset;

R = eye(size(J_DCT,2));
%%
lambda = 5e-4; %1e-4 with noise, 1e-6 or less without noise 
dctCoeff = (J_DCT'*J_DCT + lambda.^2*R)\(J_DCT'*dva);

%% inverse DCT
reconstructed_elem = DCT_subset*dctCoeff;

%% assigning element values
imgr.fwd_model = fmdl2;
imgr.fwd_model.stimulation = imdl.fwd_model.stimulation; 
imgr.elem_data = reconstructed_elem;


figure(5)
subplot(1,2,1)
show_fem(imgr);


subplot(1,2,2)
show_3d_slices(imgr, [1,1.9], [0.2],[0.5]);
view(-14,13); axis tight; axis equal;


