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
% Script overview: This is a demonstration on how to use DCT-based EIT
% algorithm. EIDORS tool bos should be set up before using this code.

load('dct_demonstration.mat')

% calculating the DCT subset with frequencies combination, which is the size to
% the prior
[dct_p, dct_q] = create_dct_subset(size(prior_l,1),size(prior_l,2));

% calculating the T mapping
specMtxCol = cal_mapping(dct_p, dct_q, fmdl, prior_l);

% check if the columns of specMtxCol are linearly independent
A = rref(specMtxCol);
s = sum(diag(A));
n = size(A,2);

if n>s
    disp('vectors are linearly dependent')
else
    disp('vectors are linearly independent')
end

% solve the inverse problem
imgRec = inv_solve_DCT(imdl, deltaVolt, specMtxCol, lambda);

% display the image
show_fem(imgRec);