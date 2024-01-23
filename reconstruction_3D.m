% from 
% https://www.sce.carleton.ca/faculty/adler/eidors/tutorial/EIDORS_basics/basic_3d.shtml

close all

% Basic 3d model $Id: basic_3d_01.m 2161 2010-04-04 20:33:46Z aadler $

fmdl= ng_mk_cyl_models(3,[15,1,1.5,2],[0.1,0,0.05]); 
show_fem(fmdl);

imdl = mk_common_model('a2c2',8); % Will replace most fields
imdl.fwd_model = fmdl;
imdl.fwd_model.stimulation = mk_stim_patterns(45,1,[0,3],[0,1],{},1);
img1 = mk_image(imdl);

figure(1)
subplot(1,2,1)
show_fem(img1); axis tight;

% print_convert('basic_3d_01a.png','-density 60')

%%
% Basic 3d model $Id: basic_3d_02.m 3790 2013-04-04 15:41:27Z aadler $

% Add a circular object at 0.2, 0.5
% Calculate element membership in object
select_fcn = inline('(x-0.2).^2 + (y-0.5).^2 + (z-2).^2 < 0.3^2','x','y','z');
memb_frac = elem_select( img1.fwd_model, select_fcn);
img2 = mk_image(img1, 1 + memb_frac );

img2.calc_colours.cb_shrink_move = [0.3,0.6,0.02];

subplot(1,2,2)
show_fem(img2,1); axis tight;

% print_convert('basic_3d_02a.png','-density 60');

%%
% Basic 3d model $Id: basic_3d_03.m 2161 2010-04-04 20:33:46Z aadler $


% Show 3D object as slices
img2.calc_colours.greylev = -0.05;
img2.calc_colours.npoints = 128;

figure(2)
subplot(1,2,1)
show_3d_slices(img2, [0.5,1.5,1.8,2.1]);
view(-14,13); axis tight; axis equal; zlim([0,3]);

% print_convert('basic_3d_03a.png','-density 60')
subplot(1,2,2)
show_3d_slices(img2, [1,1.9], [0.5],[0.5]);
view(-14,13); axis tight; axis equal;

% print_convert('basic_3d_03b.png','-density 60')

%%
% Basic 3d model $Id: basic_3d_04.m 2161 2010-04-04 20:33:46Z aadler $

% Simulate Voltages and plot them
vh= fwd_solve(img1);
vi= fwd_solve(img2);

figure(3)
plot([vh.meas, vi.meas]);
axis tight
% print_convert('basic_3d_04a.png','-density 60',0.4);

%%
% Reconstruction Model $Id: basic_3d_05.m 3126 2012-06-08 16:17:56Z bgrychtol $
J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
iRtR = inv(prior_noser( imdl ));
hp = 0.17;
iRN = hp^2 * speye(size(J,1));
RM = iRtR*J'/(J*iRtR*J' + iRN);
imdl.solve = @solve_use_matrix; 
imdl.solve_use_matrix.RM  = RM;

%%

% Reconstruct Model $Id: basic_3d_06.m 2161 2010-04-04 20:33:46Z aadler $
imgr = inv_solve(imdl, vh, vi);

imgr.calc_colours.ref_level = 0; % difference imaging
imgr.calc_colours.greylev = -0.05;

figure(4)
subplot(1,2,1)
show_fem(imgr);
% print_convert('basic_3d_06a.png','-density 60');

subplot(1,2,2)
show_3d_slices(imgr, [1,1.9], [0.5],[0.5]);
view(-14,13); axis tight; axis equal;
% print_convert('basic_3d_06b.png','-density 60');