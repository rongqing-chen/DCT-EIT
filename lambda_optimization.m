function [lambdas, logging] = lambda_optimization(voltage_homo, voltage_inho, target_noise_figure, eidors_img, subset, parameters)
% [optimal_lambda, non_single_solution_flag] = lambda_optimization(voltage_homo, 
%   voltage_inho, target_noise_figure, eidors_img, subset, parameters)
%   finds the optimal_lambda that matches the target_noise_figure. 
%   It throws a warning in case it finds multiple lambdas and updates 
%   problem_flag to true.
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
problem_flag.no_solution = false;
problem_flag.multiple_solutions = false;

%% lets keep a log
logging = struct('lambdas', [], 'decade', [], 'repetitions', []);

%% Hyperparameter range
% To improve numerical stability the evaluations points are clustered at
% the edges and there is some averaging.

h1 = parameters.lambda_low;
h2 = parameters.lambda_high;
n_points = parameters.lambda_points;
n_averages = parameters.n_average;

n_decades = max(round(log10(h2)-log10(h1)),1);

[hypervalues_1_dec, x_lin] = log_lobato_points(h1, h1*10, n_points);

noise_figures = zeros(length(hypervalues_1_dec), n_averages, n_decades);

%% Calculating the noise figure in every decade
J = calc_jacobian(eidors_img);
vi_noise = voltage_inho;

n_points_polyval = 5*n_points;
all_hypervalues = zeros(length(hypervalues_1_dec), n_decades);
all_polyval = zeros(n_points_polyval, n_decades); 
all_x_logscale = zeros(n_points_polyval, n_decades);

for idec = 1:n_decades
    hypervalues = hypervalues_1_dec*(10.^(idec-1));
    logging(idec).decade = [min(hypervalues), max(hypervalues)];
    all_hypervalues(:,idec) = hypervalues;
    
    optimal_lambda = [];
    curve_to_fit = -inf;
    max_repetitions = parameters.max_repetitions+1;
    repeat_for_sure = 0;
    while repeat_for_sure <= max_repetitions && ...
            isempty(optimal_lambda) && ...
            any(curve_to_fit <= target_noise_figure/5)
        repeat_for_sure = repeat_for_sure +1;

        for ival = 1:length(hypervalues_1_dec)
            for iavg = 1:n_averages
                noise = std(voltage_inho.meas - voltage_homo.meas)*noise_ampl;
                vi_noise.meas = voltage_inho.meas + noise * randn(size(voltage_inho.meas));
                
                % calc voltage difference
                delta_volt = calc_difference_data(voltage_homo, voltage_inho, fwd_model);
                delta_volt_noise = calc_difference_data(voltage_homo, vi_noise, fwd_model);
                
                system_noise = delta_volt - delta_volt_noise;
                
                J_subset = J* subset;
            
                lambda = hypervalues(ival);
                
                [noise_figure, ~] = calc_noise_figure(delta_volt, system_noise, J_subset, subset, lambda);
           
                noise_figures(ival,iavg,idec) =  noise_figure;        
            end
        end


    %% Find the optimal hyperparameter
    % Shift the noise figure curve by the target_noise_figure and fit with a
    % polynomial and then look for the possibly only root in the fitting 
    % interval. 
    % Some tricks: i) cluster fitting points at the edges, ii) use
    % linear (non-physcial) space for x
    
        curve_to_fit = mean(squeeze(noise_figures(:,:,idec)),2) - target_noise_figure;
        
        P = polyfit(x_lin, curve_to_fit, parameters.polynomial_order);
        h1 = log10(min(hypervalues));
        h2 = log10(max(hypervalues));
        all_x_logscale(:,idec) = logspace(h1, h2, n_points_polyval);
        x_lin_ext = linspace(-1, 1, n_points_polyval);
        all_polyval(:,idec) = polyval(P, x_lin_ext);
        
        % find the right root (hopefully only one) in the linear (non-physical) space
        lin_roots = roots(P);
        l = find((abs(imag(lin_roots))<1e-15) & (abs(real(lin_roots)) < 1));
        right_root = real(lin_roots(l));
        
        % shift the root to the physical space 
        a = h1;
        b = h2;
        
        optimal_lambda = 10.^(0.5*(a+b) + 0.5*(b-a)*right_root);
        logging(idec).lambdas = optimal_lambda;            
        logging(idec).repetitions = repeat_for_sure;

%% Plot if you need
        if parameters.make_plot
            if idec == 1 && repeat_for_sure == 1
                figure(1000)
                clf
            end
        
            hold on
            plot(all_hypervalues(:,idec), mean(squeeze(noise_figures(:,:,idec)),2), 'b')
            plot(all_x_logscale(:,idec), all_polyval(:,idec) + target_noise_figure, 'r')
            plot(logging(idec).lambdas, target_noise_figure*ones(size(logging(idec).lambdas)), 'ko')
            
            line(xlim,[target_noise_figure,target_noise_figure], 'color', 'k')

            set(gca, 'xscale', 'log')
            xlabel('Lambda')
            ylabel('Noise Figure')
        end
            
    end

    if length(optimal_lambda) >= 1 && parameters.break_fast % finished
        break
    end

end    

%% check
% lambdas = cell2mat({logging.lambdas});
% lambdas(isnan(lambdas)) = [];

lambdas = get_lambdas_from_struct(logging);

if isempty(lambdas)
    warning('!!No solution found!!')
    problem_flag.no_solution = true;
elseif length(lambdas) > 1
    warning('!!Multiple solutions found!!')
    problem_flag.multiple_solutions = true;
end

end


%%
function values = get_lambdas_from_struct(logging)
values = [];
idx = 0;
for ii = 1:length(logging)
    if ~isempty(logging(ii).lambdas)
        idx = idx +1;
        val = logging(ii).lambdas;
        values(idx:idx+length(val)-1) = val;
    end
end

end


function [updated_parameters, default_parameters] = my_default_parameters(parameters)
    % Define a structure with default values
    default_parameters = struct('noise_amplitude', 0.01, ...
        'lambda_low', 1e-3,...
        'lambda_high', 1e-1, ...
        'lambda_points', 7, ...
        'polynomial_order', 3, ...
        'n_average', 1, ...
        'max_repetitions', 5, ...
        'break_fast', true, ...
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




