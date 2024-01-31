function [optimal_lambda, non_single_solution_flag] = lambda_optimization(voltage_homo, voltage_inho, target_noise_figure, eidors_img, subset, parameters)
% [optimal_lambda, non_single_solution_flag] = lambda_optimization(voltage_homo, 
%   voltage_inho, target_noise_figure, eidors_img, subset, parameters)
%   finds the optimal_lambda that matches the target_noise_figure. 
%   It throws a warning in case it finds multiple lambdas and updates 
%   non_single_solution_flag to true.
%   parameters have some default values.

arguments
    voltage_homo
    voltage_inho
    target_noise_figure
    eidors_img
    subset
    parameters = [];
end

% update the parameters with the defualt if needed
parameters = my_default_parameters(parameters);

noise_ampl = parameters.noise_amplitude;

fwd_model = eidors_img.fwd_model;

% hopefully everything goes fine
non_single_solution_flag = false;


%% Hyperparameter range
% To improve numerical stability the evaluations points are clustered at
% the edges and there is some averaging.

h1 = parameters.lambda_low;
h2 = parameters.lambda_high;
n_points = parameters.lambda_points;
n_averages = parameters.n_average;
[hypervalues, x_lin] = log_lobato_points(h1, h2, n_points);

noise_figures = zeros(length(hypervalues), n_averages);


%% Calculating the noise figure
J = calc_jacobian(eidors_img);
vi_noise = voltage_inho;

for jj = 1:length(hypervalues)
    for ii = 1:n_averages
        noise = std(voltage_inho.meas - voltage_homo.meas)*noise_ampl;
        vi_noise.meas = voltage_inho.meas + noise * randn(size(voltage_inho.meas));
        
        % calc voltage difference
        delta_volt = calc_difference_data(voltage_homo, voltage_inho, fwd_model);
        delta_volt_noise = calc_difference_data(voltage_homo, vi_noise, fwd_model);
        
        system_noise = delta_volt - delta_volt_noise;
        
        J_subset = J* subset;
    
        lambda = hypervalues(jj);
        
        [noise_figure, ~] = calc_noise_figure(delta_volt, system_noise, J_subset, subset, lambda);
   
        noise_figures(jj,ii) =  noise_figure;        
    end
end


%% Find the optimal hyperparameter
% Shift the noise figure curve by the target_noise_figure and fit with a
% polynomial and then look for the possibly only root in the fitting 
% interval. 
% Some tricks: i) cluster fitting points at the edges, ii) use
% linear (non-physcial) space for x

curve = mean(noise_figures,2) - target_noise_figure;

P = polyfit(x_lin, curve, 7);

x_logscale = logspace(log10(h1), log10(h2), 5*n_points);
x_lin_ext = linspace(-1, 1, 5*length(hypervalues));
y = polyval(P, x_lin_ext);

% find the right root (hopefully only one) in the linear (non-physical) space
lin_roots = roots(P);
l = find((abs(imag(lin_roots))<1e-15) & (abs(real(lin_roots)) < 1));
right_root = real(lin_roots(l));

% shift the root to the physical space 
a = log10(h1);
b = log10(h2);

optimal_lambda = 10.^(0.5*(a+b) + 0.5*(b-a)*right_root);

if length(optimal_lambda) > 1 % ops
    warning('!!More than a possible lambda!!')
    non_single_solution_flag = true;
end


%% Plot if you need
if parameters.make_plot
    figure(1000)
    clf
    hold on
    plot(hypervalues, curve + target_noise_figure)
    plot(x_logscale, y + target_noise_figure)
    plot(optimal_lambda, target_noise_figure*ones(size(optimal_lambda)), 'o')
    
    set(gca, 'xscale', 'log')
    xlabel('Lambda')
    ylabel('Noise Figure')
end

end

%%
function [updated_parameters, default_parameters] = my_default_parameters(parameters)
    % Define a structure with default values
    default_parameters = struct('noise_amplitude', 0.01, ...
        'lambda_low', 1e-3,...
        'lambda_high', 1e-1, ...
        'lambda_points', 20, ...
        'n_average', 2, ...
        'make_plot', false);
    
    updated_parameters = default_parameters;
    
    if isempty(parameters)
        return
    end

    % Update updated_parameters with values from parameters, if provided
    param_names = fieldnames(default_parameters);
    for i = 1:length(param_names)
        param_name = param_names{i};
        if isfield(parameters, param_name)
            updated_parameters.(param_name) = parameters.(param_name);
        end
    end

end




