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
% Script overview: 
% This is a function to create DCT basic cosine subsets
% It does NOT include any prior information inside
% The function returns a DCT subset which is used for FEM clustering
% X_Pix by Y_Pix is the size (pixel by pixel) of given image

function [dct_p, dct_q] = create_dct_subset(X_Pix, Y_Pix)

x_range = 0:X_Pix - 1;
y_range = 0:Y_Pix - 1;

[pp,qq] = meshgrid(x_range, y_range);

% MD DCT II
dct_p = sqrt(2 / X_Pix) * cos(pi * (2*pp + 1) .* qq / (2 * X_Pix));
dct_p(1,:) = dct_p(1,:) / sqrt(2);

dct_q = sqrt(2 / Y_Pix) * cos(pi * (2*pp + 1) .* qq / (2 * Y_Pix));
dct_q(1,:) = dct_q(1,:) / sqrt(2);

end