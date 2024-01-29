function [subset] = make_subset_3D(elem_centers, coefficients_matrix, subset_makers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% subset = zeros(length(elem_centers), length(coefficients_matrix));
single_dim_subsets = zeros(length(elem_centers), length(coefficients_matrix), 3);
for idm = 1:3
    single_dim_subsets(:,:,idm) = subset_makers{idm}(elem_centers(:,idm),coefficients_matrix(:,idm));
end

subset = prod(single_dim_subsets,3);


% col = 0;
% for ix = 1:size(single_dim_subsets{1,1},2)
%     for iy = 1:size(single_dim_subsets{2,1},2)
%         for iz = 1:size(single_dim_subsets{3,1},2)
%             col = col+1;
%             subset(:,col) = single_dim_subsets{1,1}(:,ix).*single_dim_subsets{2,1}(:,iy).*single_dim_subsets{3,1}(:,iz);
%         end
%     end
% end





% for col = 1:length(coefficients_matrix)
%     a_jk = coefficients_matrix(col,:);
%     dims = zeros(length(elem_centers),3);
%     for idm = 1:3
%         dims(:,idm) = subset_makers{1}(elem_centers(:,idm),a_jk(idm));
%     end
% 
%     subset(:,col) = prod(dims,2);
% 
% end


end