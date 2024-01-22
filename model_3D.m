
clear
home

% from 
% https://www.sce.carleton.ca/faculty/adler/eidors/tutorial/EIDORS_basics/forward_solvers_3d.shtml

nelec= 8; nrings= 2;
ring_vert_pos = [0.2, 0.5]; 
fmdl= ng_mk_cyl_models([1,0.3,0.05],[nelec,ring_vert_pos],[0.1,0.05,0.02]);

stim = mk_stim_patterns(nelec,nrings,[0,1],[0,1],{'meas_current'},1);
fmdl.stimulation = stim;

conduct = 1;
img = mk_image( fmdl, conduct ); 

show_fem(img);
% print_convert forward_solver_3d_01a.png '-density 75'

%%
extra={'ball','solid ball = sphere(0,0.2,0.5;0.1);'};
fmdl= ng_mk_cyl_models([1,0.3,0.05],[nelec,ring_vert_pos],[0.1,0.05,0.02], extra);
fmdl.stimulation = stim;

img= mk_image(fmdl, conduct);
img.elem_data(fmdl.mat_idx{2}) = 0.1;

show_fem(img);
% print_convert forward_solver_3d_02a.png '-density 75'

%%
img.fwd_solve.get_all_meas = 1;
img.elem_data(fmdl.mat_idx{2}) = conduct; % Homogenous
vh = fwd_solve(img);

img.elem_data(fmdl.mat_idx{2}) = 0.95*conduct; %Non-conductive inclusion
vi = fwd_solve(img);

plot([vh.meas, 100*(vh.meas-vi.meas)])
axis tight
% print_convert forward_solver_3d_03a.png '-density 77'

%%
img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);

show_slices(img_v,[0.2;0.3;0.4]*[inf,inf,1])
% print_convert forward_solver_3d_04a.png '-density 75'

img_v.node_data = vh.volt(:,1) - vi.volt(:,1);
show_slices(img_v,[0.2;0.3;0.4]*[inf,inf,1])
% print_convert forward_solver_3d_04b.png '-density 75'

%%
img.fwd_model.electrode([2,13]) = img.fwd_model.electrode([13,2]); % flip electrodes

img.elem_data(fmdl.mat_idx{2}) = conduct; % Homogenous
vh = fwd_solve(img);

img.elem_data(fmdl.mat_idx{2}) = 0.95*conduct; %Non-conductive inclusion
vi = fwd_solve(img);

img_v.node_data = vh.volt(:,1);

show_slices(img_v,[0.25;0.35;0.45]*[inf,inf,1])
% print_convert forward_solver_3d_05a.png '-density 75'

img_v.node_data = vh.volt(:,1) - vi.volt(:,1);
show_slices(img_v,[0.25;0.35;0.45]*[inf,inf,1])
% print_convert forward_solver_3d_05b.png '-density 75'