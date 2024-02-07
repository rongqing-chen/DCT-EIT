%% plot singular values for J and J_subset

s = svd(J);
s_subset = svd(J_subset);

figure(1)
clf
hold on
plot(s(1:end/2)./max(s))
plot(s_subset(1:end/2)./max(s_subset))
hold off

set(gca, 'xscale', 'log', 'yscale', 'log')