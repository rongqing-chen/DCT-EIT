
% init_eidors()
clear
% 
load('dct_demonstration.mat')

% prior_l = ones(size(prior_l));

%%
% calculating the DCT subset with frequencies combination, which is the size to
% the prior
[dct_p, dct_q] = create_dct_subset(size(prior_l,1),size(prior_l,2));


%% scale model

% scale and shift model
magic_values(1,:) = [8, 38.72]; % orig_mins
magic_values(2,:) = [248, 217.28]; % orig_maxs;
magic_values(3,:) = [120, 256/2]; % stretch and shift
[fmdl_stretch, new_centers] = scale_model_dimension(fmdl, magic_values);

%% my dct
% number of coefficients
M = 16;
N = 16;

% % coefficients ordered in natural way
% [MM, NN] = ndgrid(0:M-1, 0:N-1);
% coefficients_matrix = [MM(:), NN(:)];

% order coefficients in zig zag way
coefficients_matrix = order_coeffs_tensor_product(0:M-1, 0:N-1);

S_my_dct = make_DCT_subset(new_centers, coefficients_matrix);


%% polynomial subset
% number of coefficients
M = 16;
N = M;
coefficients_matrix = order_coeffs_tensor_product(0:M-1, 0:N-1);
[norm_centers] = (2*recenter_model(fmdl))-1;
S_lege = make_legendre_subset(norm_centers, coefficients_matrix);
S_cheb = make_chebyshev_subset(norm_centers, coefficients_matrix);



%% mask

unstruct_maks = make_unstructured_mask(fmdl_stretch, prior_l);

% apply mask on subset
masked_my_dct = unstruct_maks.*S_my_dct;
masked_lege = unstruct_maks.*S_lege;
masked_cheb = unstruct_maks.*S_cheb;


% check OK
imgRec_my_dct = mk_image(imdl,1);
imgRec_my_dct.elem_data = unstruct_maks;

figure(1)
clf
subplot(2,1,1)
show_fem(imgRec_my_dct)

subplot(2,1,2)
imagesc(flipud(prior_l))
% axis equal


%%
% calculating the T mapping
spec_Mtx_col = cal_mapping(dct_p, dct_q, fmdl, prior_l);


%%
% solve the inverse problem
imgRec = inv_solve_DCT(imdl, deltaVolt, spec_Mtx_col, lambda);


%% my version
imgRec_my_dct = inv_solve_DCT(imdl, deltaVolt, masked_my_dct, lambda);
imgRec_lege = inv_solve_DCT(imdl, deltaVolt, masked_lege, lambda);
imgRec_cheb = inv_solve_DCT(imdl, deltaVolt, masked_cheb, lambda);


%%
% display the image
figure(2)
clf
subplot(2,2,1)
show_fem(imgRec);

subplot(2,2,3)
show_fem(imgRec_my_dct);

subplot(2,2,2)
show_fem(imgRec_lege);

subplot(2,2,4)
show_fem(imgRec_cheb);

%% compare singular values
% I am not sure about this part
imgRec_ = mk_image(imdl,1);

% calc the Jacobian
J = calc_jacobian(imgRec_);

S_J = svd(J); 
S_dct = svd(J*spec_Mtx_col);
S_my_dct = svd(J*masked_my_dct); 
S_lege = svd(J*masked_lege); 
S_cheb = svd(J*masked_cheb); 

figure(3)
clf
hold on
loglog(S_J./S_J(1))
loglog(S_dct./S_dct(1))
loglog(S_my_dct./S_my_dct(1))
loglog(S_lege./S_lege(1))
loglog(S_cheb./S_cheb(1))

set(gca, 'xscale', 'lin',  'yscale', 'log')

legend('J', 'dct', 'my dct', 'lege', 'cheb')


%%
cond(J)
cond(J.*unstruct_maks')
cond(J*spec_Mtx_col)
cond(J*masked_my_dct)
cond(J*masked_lege)
cond(J*masked_cheb)