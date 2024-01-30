% %% This is a demonstration of the NF calculation 
% 
% extra={'ball','solid ball = sphere(0,0.2,1.2;0.5);'};
% fmdl = ng_mk_cyl_models([3,1,0.2],[16,1.5],[0.1,0,0.05], extra); % 2994 nodes
% imdl = mk_common_model('a2c2',8); % Will replace most fields
% imdl.fwd_model = fmdl;
% stim_pattern = mk_stim_patterns(16,1,[0,3],[0,1],{},1);
% imdl.fwd_model.stimulation = stim_pattern;
% 
% img = mk_image(imdl);
% vh = fwd_solve(img);
% 
% img_extra = img;
% img_extra.elem_data(fmdl.mat_idx{2}) = 0.1;
% 
% vi = fwd_solve(img_extra);
% %%
% 
% add 25% noise to the data
vi_noise = vi;
noise_ampl = 0.1;
noise = std(vi.meas - vh.meas)*noise_ampl;
vi_noise.meas = vi.meas + noise_ampl * randn(size(vi.meas));

% calc voltage difference
delta_volt = calc_difference_data(vh, vi, imdl.fwd_model);
delta_volt_noise = calc_difference_data(vh, vi_noise, imdl.fwd_model);

system_noise = delta_volt - delta_volt_noise;
% 
% %%
% img_rec = mk_image(imdl_rec);
% 
% 
% %% subset jacobian
% J = calc_jacobian(img_rec);

%% This part you can define your own range of your hyperparameter
hypervalues = logspace(-3,-1,30);

imgRec.calc_colours.ref_level =  0;

numb_averages = 10;
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

figure(10)
semilogx(hypervalues, mean(NF_DCT,2))

%% how to use the function in EIDORS
% imdl.hyperparameter.func = @choose_noise_figure;
% imdl.hyperparameter.noise_figure = 0.5;
% % very tricky here, you need to assign a vector of element numbers of
% % contrast in centre as a simulation target
% imdl.hyperparameter.tgt_elems = [2037,2038];
% choose_noise_figure(imdl); %calc_noise_figure is called by choose_noise_figure


