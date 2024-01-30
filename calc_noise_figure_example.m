
noise_ampl = 0.01;


%% This part you can define your own range of your hyperparameter
% hypervalues = logspace(-3,-1,20);
h1 = 1e-3;
h2 = 1e-1;
n_points = 20;
[hypervalues, x_lin] = log_lobato_points(h1, h2, n_points);


imgRec.calc_colours.ref_level =  0;

target_noise_figure = 0.5;

numb_averages = 2;
NF_DCT = zeros(length(hypervalues), numb_averages);


%% calculating the noise figure
for j = 1:length(hypervalues)
    for ii = 1:numb_averages
        noise = std(vi.meas - vh.meas)*noise_ampl;
        vi_noise.meas = vi.meas + noise * randn(size(vi.meas));
        
        % calc voltage difference
        delta_volt = calc_difference_data(vh, vi, imdl.fwd_model);
        delta_volt_noise = calc_difference_data(vh, vi_noise, imdl.fwd_model);
        
        system_noise = delta_volt - delta_volt_noise;
    
    
        J_DCT = J* DCT_subset;
    
        lambda = hypervalues(j);
        
        [noise_figure, reconstructed_elem] = calc_noise_figure(delta_volt, system_noise, J_DCT, DCT_subset, lambda);

    
        NF_DCT(j,ii) =  noise_figure;
        
    end
end

%%

curve = mean(NF_DCT,2) - target_noise_figure;

P = polyfit(x_lin, curve, 7);

x_logscale = logspace(log10(h1), log10(h2), 5*n_points);
x_lin_ext = linspace(-1, 1, 5*length(hypervalues));
y = polyval(P, x_lin_ext);

% 
lin_roots = roots(P);
l = find((abs(imag(lin_roots))<1e-15) & (abs(real(lin_roots)) < 1));
right_root = real(lin_roots(l));

a = log10(h1);
b = log10(h2);

lambda_0 = 10.^(0.5*(a+b) + 0.5*(b-a)*right_root);

if length(lambda_0) > 1
    disp('!!more than a solution!!')
else
    disp('OK')
end

%
figure(10)
clf
hold on
plot(hypervalues, curve + target_noise_figure)
plot(x_logscale, y + target_noise_figure)
plot(lambda_0, target_noise_figure*ones(size(lambda_0)), 'o')

set(gca, 'xscale', 'log')





