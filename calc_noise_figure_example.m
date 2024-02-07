
parameters.noise_amplitude = 0.01;
parameters.lambda_low = 1e-3;
parameters.lambda_high = 1e-1;
parameters.lambda_points = 20;
parameters.n_average = 2;
parameters.make_plot = true;

target_noise_figure = 0.5;
eidors_img = img_rec;
subset = DCT_subset;

unsuccessful_times = 0;
lambdas = [];
runs = 0;
repetitions = 1000;
for ii = 1:repetitions
    [lambda_0, non_single_solution_flag] = lambda_optimization(vh, vi, target_noise_figure, eidors_img, subset, parameters);
    if non_single_solution_flag
        unsuccessful_times = unsuccessful_times+1;
    else
        runs = runs +1;
        lambdas(runs) = lambda_0;
    end
end

%%
fprintf('The procedure was unsuccessful %.1f percent of the times\n', unsuccessful_times/repetitions*100)
mean_lambdas = mean(lambdas);
std_lambdas = std(lambdas);

fprintf('Optimal lambda = %.3e +- %.3e \n', mean_lambdas, 2*std_lambdas)

y = linspace(mean_lambdas-3*std_lambdas, mean_lambdas+3*std_lambdas, 200);
f = exp(-(y-mean_lambdas).^2./(2*std_lambdas^2))./(std_lambdas*sqrt(2*pi));

figure(10)
clf
histogram(lambdas, max(round(repetitions/20),5), 'Normalization', 'pdf')
hold on
plot(y,f,'LineWidth',1.5)

xlabel('Optimal Lambda')
ylabel('Probability Density')

% The procedure was unsuccessful 1.2 percent of the times
% Optimal lambda = 5.497e-03 +- 1.081e-03 
% lambda distribution is asymmetric
