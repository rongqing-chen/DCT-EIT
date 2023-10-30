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
% Script overview: This function reconstruct the an image of conductivity
% change using the oundary measurement data. An inverse model different to
% the forward model should be used to avoid inverse crime. Lambda should be
% chosen using proper method, such as NF.
% 

function imgRec = inv_solve_DCT(imdl, deltaVolt, specMtxCol, lambda)

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

end
