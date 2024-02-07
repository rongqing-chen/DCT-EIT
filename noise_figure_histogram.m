%% Show that the noise_figures are normally distributed
close all
hypervalues = logspace(-3,-1,6);

numb_averages = 500;
noise_figures = zeros(length(hypervalues), numb_averages);

vi_noise = vi;

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

    
        noise_figures(j,ii) =  noise_figure;
        
    end
end

%%
figure(1)
clf
tiledlayout('flow')
for j = 1:length(hypervalues)
    nexttile
    histogram(noise_figures(j,:),30,'Normalization','pdf')
    mu = mean(noise_figures(j,:));
    sigma = std(noise_figures(j,:));
    y = linspace(mu-3*sigma,mu+3*sigma, 200);
    f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    
    hold on
    plot(y,f,'LineWidth',1.5)
end