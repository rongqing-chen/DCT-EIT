

ll = [10,125,200]; 

for iii = 1:length(ll)
    figure(iii)
    clf
    hold on
    plot(spec_Mtx_col(:,ll(iii))/max(abs(spec_Mtx_col(:,ll(iii)))))
    plot(values(:,ll(iii))/max(abs(values(:,ll(iii)))))
    hold off
end
