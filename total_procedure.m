
% init_eidors()
clear
% 
load('dct_demonstration.mat')

% prior_l = ones(size(prior_l));

%%
% calculating the DCT subset with frequencies combination, which is the size to
% the prior
% [dct_p, dct_q] = create_dct_subset(size(prior_l,1),size(prior_l,2));

[X_Pix, Y_Pix] = size(prior_l);
x_range = 0:X_Pix - 1;
y_range = 0:Y_Pix - 1;

[pp,qq] = meshgrid(x_range, y_range);

dct_p = sqrt(2 / X_Pix) * cos(pi * (2*pp + 1) .* qq / (2 * X_Pix));
dct_p(1,:) = dct_p(1,:) / sqrt(2);

dct_q = sqrt(2 / Y_Pix) * cos(pi * (2*pp + 1) .* qq / (2 * Y_Pix));
dct_q(1,:) = dct_q(1,:) / sqrt(2);


%% scale model

% scale and shift model
magic_values(1,:) = [8, 38.72]; % orig_mins
magic_values(2,:) = [248, 217.28]; % orig_maxs;
magic_values(3,:) = [120 , 256/2]; % stretch and shift
[fmdl_stretch, new_centers] = scale_model_dimension(fmdl, magic_values);

%% my dct
% number of coefficients
M = 16;
N = 16;


[S, ordered_coefficients] = make_DCT_subset(new_centers, M, N);


%% mask
pts = interp_mesh(fmdl_stretch, 0);

x = 0:255;
y = x';

mask_interp = griddedInterpolant({x,y}, prior_l, 'nearest');

% dimensions need inversion
unstruct_maks = mask_interp([pts(:,2), pts(:,1)]);

%
masked_values = unstruct_maks.*S;


% check OK
imgRec = mk_image(imdl,1);
imgRec.elem_data = unstruct_maks;

figure(1)
clf
subplot(2,1,1)
show_fem(imgRec)

subplot(2,1,2)
imagesc(flipud(prior_l))
% axis equal

%%
% calculating the DCT subset with frequencies combination, which is the size to
% the prior
% [dct_p, dct_q] = create_dct_subset(size(prior_l,1),size(prior_l,2));



% preparing the prior information
prior = flipud(prior_l);

% recenter and stretch the model to 256x256
fmdl_stretch = fmdl;
fmdl_stretch.nodes = fmdl.nodes * 120 + 256/2;
pts = interp_mesh(fmdl_stretch, 0); % center of elements
round_pts = round(pts);

%
% img = mk_image(fmdl,1);
% image_values = prior_l(round_pts);
% img.elem_data = 1 - 0.8*image_values;

temp_idx = sub2ind(size(prior), round_pts(:,2), round_pts(:,1));

% get the dimension at the x-axis and y-axis
dimX = sqrt(length(dct_p));
dimY = sqrt(length(dct_q));

spec_Mtx = zeros(size(prior_l,1),size(prior_l,2),dimX*dimY);
% spec_Mtx_col = zeros(size(temp_idx,1),dimX*dimY);

% create subset
index = 0;
for i=1:dimX
    for j = 1:dimY
        index = index + 1;
        temp_Mtx = dct_p(i,:)'*dct_q(j,:);
        spec_Mtx(:,:,index) = temp_Mtx;
    end
end

% apply mask
masked_spec_Mtx = reshape(spec_Mtx.*prior_l,[],dimX*dimY);
spec_Mtx_col = masked_spec_Mtx(temp_idx, :);


%%
specMtxCol = spec_Mtx_col;



%%
% solve the inverse problem
% function imgRec = inv_solve_DCT(imdl, deltaVolt, specMtxCol, lambda)

% creating an inverse model for EIT

imgRec = mk_image(imdl,1);

% calc the Jacobian
J = calc_jacobian(imgRec);

%% solving the change of DCT coefficients
J_DCT = J * specMtxCol;
[~, p] = size(J_DCT);
R = eye(p);

dctCoeff = (J_DCT'*J_DCT + lambda.^2*R)\(J_DCT'*deltaVolt);

%% inverse DCT

% recCond = zeros(size(imgRec.elem_data));

% for i=1:numel(dctCoeff)
%     recCond = recCond + specMtxCol(:,i) .* dctCoeff(i);
% end

recCond = specMtxCol*dctCoeff;

%% my version
J_my_DCT = J * masked_values;
R = eye(size(J_my_DCT,2));

my_dctCoeff = (J_my_DCT'*J_my_DCT + lambda.^2*R)\(J_my_DCT'*deltaVolt);

%% my inverse DCT
my_recCond = masked_values*my_dctCoeff;


%% assigning element values

imgRec.elem_data = recCond;
my_imgRec = imgRec;

my_imgRec.elem_data = my_recCond;

%%
% display the image
figure(2)
clf
subplot(2,1,1)
show_fem(imgRec);

subplot(2,1,2)
show_fem(my_imgRec);

%% compare singular values
S_J = svd(J);
S_dct = svd(J*spec_Mtx_col);
S_val = svd(J*S); 

figure(3)
clf
hold on
loglog(S_J./S_J(1))
loglog(S_dct./S_dct(1))
loglog(S_val./S_val(1))

set(gca, 'xscale', 'lin',  'yscale', 'log')

legend('J', 'dct', 'my dct')