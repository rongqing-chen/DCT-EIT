
% init_eidors()

% 
load('dct_demonstration.mat')

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


%%
% calculating the DCT subset with frequencies combination, which is the size to
% the prior
% [dct_p, dct_q] = create_dct_subset(size(prior_l,1),size(prior_l,2));

% preparing the prior information
prior = flipud(prior_l);

% recenter and stretch the model to 256x256
fmdlStretch = fmdl;
fmdlStretch.nodes = fmdl.nodes * 120 + 256/2;
pts = interp_mesh(fmdlStretch, 0); % center of elements
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
masked_spec_Mtx = reshape(spec_Mtx.*prior_l,[],index);
spec_Mtx_col = masked_spec_Mtx(temp_idx, :);


%%
specMtxCol = spec_Mtx_col;
% check if the columns of specMtxCol are linearly independent
A = rref(specMtxCol); % not sure this is correct?!?
s = sum(diag(A));
n = size(A,2);

if n>s
    disp('vectors are linearly dependent')
else
    disp('vectors are linearly independent')
end


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

recCond = zeros(size(imgRec.elem_data));

for i=1:numel(dctCoeff)
    recCond = recCond + specMtxCol(:,i) .* dctCoeff(i);
end

%% assigning element values

imgRec.elem_data = recCond;


%%
% display the image
figure(1)
show_fem(imgRec);