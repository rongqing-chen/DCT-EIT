clear
init_eidors()

extra={'ball','solid ball = sphere(0,0.2,1.2;0.5);'};
fmdl = ng_mk_cyl_models([3,1,0.2],[16,1.5],[0.1,0,0.05], extra); % 2994 nodes
imdl = mk_common_model('a2c2',8); % Will replace most fields
imdl.fwd_model = fmdl;
stim_pattern = mk_stim_patterns(16,1,[0,3],[0,1],{},1);
imdl.fwd_model.stimulation = stim_pattern;

img = mk_image(imdl);
vh = fwd_solve(img);

img_extra = img;
img_extra.elem_data(fmdl.mat_idx{2}) = 0.1;

vi = fwd_solve(img_extra);

figure(1)
clf
subplot(1,2,1)
show_fem(img);

subplot(1,2,2)
show_fem(img_extra);

figure(2)
clf
hold on
plot(vh.meas)
plot(vi.meas)
hold off

%% model for reconstruction
fmdl_rec = ng_mk_cyl_models([3,1,0],[16,1.5],[0.1,0,0.05]); % 2143 nodes
imdl_rec = mk_common_model('a2c2',8); % Will replace most fields
imdl_rec.fwd_model = fmdl_rec;
imdl_rec.fwd_model.stimulation = stim_pattern;

img_rec = mk_image(imdl_rec);


%%
J = calc_jacobian(img_rec);
% elem_centers = (recenter_model(fmdl_rec).*[pi,2,2])-[0,1,1];
elem_centers = (recenter_model(fmdl_rec).*[pi,pi,pi])-[0,0,0];

dva = calc_difference_data( vh, vi, fmdl_rec);

M = 5;
N = 5;
O = 2;

% coefficients ordered in row, by col, by depth
[MM, NN, OO] = ndgrid(0:M-1, 0:N-1, 0:O-1);
coefficients_matrix = [MM(:), NN(:), OO(:)];

subset_makers = {@make_DCT_subset, @make_DCT_subset, @make_DCT_subset};

subset = make_subset_3D(elem_centers, coefficients_matrix, subset_makers);
DCT_subset = make_DCT_subset(elem_centers, coefficients_matrix);

J_DCT = J* DCT_subset;

R = eye(size(J_DCT,2));
%%
lambda = 5e-4; %1e-4 with noise, 1e-6 or less without noise 
disp(lambda)
dctCoeff = (J_DCT'*J_DCT + lambda.^2*R)\(J_DCT'*dva);

%% inverse DCT
reconstructed_elem = DCT_subset*dctCoeff;


img_rec.elem_data = reconstructed_elem;


figure(5)
subplot(1,2,1)
show_fem(img_rec);


subplot(1,2,2)
show_3d_slices(img_rec, [1.5], [0.2],[0.5]);
view(-14,13); axis tight; axis equal;
