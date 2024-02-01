% Visualize the subsets according to their coefficients
% run simple_3D, you need coefficients_matrix, subset, new_elem_centers,
% and img_final
% to see things clearly you should have only one of the coefficients ~=1
% and all the others =1.


close all
% combinations of coefficients which has a single coefficient ~=1
l = find(sum(coefficients_matrix==1,2)==2);
idx_coeff = 1:length(coefficients_matrix(l,:));

first_fig = 100;

for idx = idx_coeff
    subset_ = subset(:,l(idx));
    coeffs = coefficients_matrix(l(idx),:);
    img_final.elem_data = subset_;
    
    figure(first_fig+idx)
    clf
    tiledlayout('flow')
    
    nexttile
    show_fem(img_final);
    
    nexttile
    show_3d_slices(img_final, [1.5], [0.2], [0.5]);
    view(-14,13); axis tight; axis equal;
    
    nexttile
    show_slices (img_final, [inf, inf, target.center(3)] )
    
    nexttile
    show_slices (img_final, [0, inf, inf] )

    for idm = 1:3
        nexttile
        plot(new_elem_centers(:,idm), subset_, '.')
        legend(sprintf('coeff.: %d', coeffs(1,idm)))
    end

end

