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

function specMtxCol = cal_mapping(dct_p, dct_q, fmdl, prior_l)

% preparing the prior information
prior = flipud(prior_l);

% calculating the interpolation points onto mdl elements

fmdlStretch = fmdl;
fmdlStretch.nodes = fmdl.nodes * 120 + 256/2;
pts = interp_mesh(fmdlStretch, 0);

img = mk_image(fmdl,1);
transferFunc = zeros(length(pts),2);


for i=1:size(pts,1)
    imageValue = prior_l(round(pts(i,2)),round(pts(i,1)));
    % mind that the relationship between attenuation and the conductivity! 
    img.elem_data(i) = 1 - 0.8*imageValue;
    transferFunc(i,:) = [i,sub2ind(size(prior),round(pts(i,2)),round(pts(i,1)))];
end

% get the dimension at the x-axis and y-axis
dimX = sqrt(length(dct_p));
dimY = sqrt(length(dct_q));

specMtx = zeros(size(prior_l,1),size(prior_l,2),dimX*dimY);
specMtxCol = zeros(size(transferFunc,1),dimX*dimY);

% setting an index for the T mapping

index = 0;

for i=1:dimX
    for j = 1:dimY
        index = index + 1;
        specMtx(:,:,index) = dct_p(i,:)'*dct_q(j,:) .* prior_l;
        tempMtx = specMtx(:,:,index);
        tempIdx = ind2sub(size(prior),transferFunc(:,2));
        specMtxCol(:,index) = tempMtx(tempIdx);
    end
end

end
