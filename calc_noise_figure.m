function [noise_figure, reconstructed_elem] = calc_noise_figure(delta_volt, system_noise, J_subset, subset, lambda)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

R = eye(size(J_subset,2));
coeffs = (J_subset'*J_subset + lambda.^2*R)\(J_subset'*delta_volt);
coeffs_noise = (J_subset'*J_subset + lambda.^2*R)\(J_subset'*system_noise);

reconstructed_elem = subset*coeffs;
reconstructed_elem_noise = subset*coeffs_noise;

SNR_output = mean(abs(reconstructed_elem)) / std(reconstructed_elem_noise);  
SNR_input = mean(abs(delta_volt)) / std(system_noise);

noise_figure = SNR_input./SNR_output;

end