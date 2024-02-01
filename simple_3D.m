clear
% init_eidors()
target.center = [0, 0.2, 1.5];
target.radius = 0.35;
extra={'ball',sprintf('solid ball = sphere(%f,%f,%f; %f);', target.center, target.radius )};
fmdl = ng_mk_cyl_models([3,1,0.2],[16,1.5],[0.1,0,0.05], extra); % 3148 nodes
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


%% subset jacobian
J = calc_jacobian(img_rec);

elem_centers = interp_mesh(img_rec.fwd_model, 0); % center of elements

M = 2;
N = 7;
O = 3;

% coefficients ordered in row, by col, by depth
[MM, NN, OO] = ndgrid(1:M, 1:N, 1:O);
coefficients_matrix = [MM(:), NN(:), OO(:)];

[new_elem_centers, subset] = complete_subset_gen(elem_centers, [M,N,O], {'cyl', 'dct', 'poly'});


%%
target_noise_figure = 0.5;
parameters.make_plot = true;
parameters.lambda_low = 1e-3;
parameters.lambda_high = 1e3;

[lambda, logging] = lambda_optimization(vh, vi, target_noise_figure, img_rec, subset, parameters);
disp(lambda)
lambda = min(lambda);

%%
J_subset = J* subset;
R = eye(size(J_subset,2));
delta_volt = calc_difference_data( vh, vi, img_rec.fwd_model);

reconstruction_coeffs = (J_subset'*J_subset + lambda.^2*R)\(J_subset'*delta_volt);

%% inverse 
reconstructed_elem = subset*reconstruction_coeffs;

img_final = mk_image(imdl_rec);
img_final.elem_data = reconstructed_elem;


figure(5)
tiledlayout('flow')

nexttile
show_fem(img_final);

nexttile
show_3d_slices(img_final, [1.5], [0.2], [0.5]);
view(-14,13); axis tight; axis equal;

nexttile
show_slices (img_final, [inf, inf, target.center(3)] )

nexttile
show_slices (img_final, [target.center(1), inf, inf] )
