% Copyright 2019 Province of British Columbia
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Author: Rongqing Chen
% Date: 13-Oct-2021
%
% Script overview: This function assembles the T mapping function in DCT
% algorithm. The DCT subset, the forward model and the prior information
% must be given.
% 

function spec_Mtx_col = cal_mapping(dct_p, dct_q, fmdl, prior_l)

% preparing the prior information
prior = flipud(prior_l);

% recenter and stretch the model to 256x256
fmdlStretch = fmdl;
fmdlStretch.nodes = fmdl.nodes * 120 + 256/2;
pts = interp_mesh(fmdlStretch, 0);
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

end
