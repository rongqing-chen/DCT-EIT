
J_my_DCT = J * values;

J_2_check = {J, J_DCT, J_my_DCT};
labels = {'J', 'DCT', 'my DCT'};

inv_fun = @pinv; % @pinv or @inv, but it gives conditioning problems

figure(1)
clf
hold on
for ii = 1:length(J_2_check)
    j = J_2_check{ii}; 
    [Q,R] = qr(j);
    s = sort(diag(inv_fun(R'*R)), 'descend');
    plot(s/s(1), 'DisplayName', labels{ii})
end

set(gca, 'yscale', 'log', 'xscale', 'log')
legend