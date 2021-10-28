L = 1e5 ; % -> the number of SNR samples
N = 128 ; % -> the RIS size
P_max_dBm = 30 ; % -> the maximum power of the PA at the RIS in dBm
P_t_dBm = 30 ; % -> transmit power at the BS in dBm
G_max_dB = 30 ; % -> maximum gain of the PA in dB
n_amp_dBm = 100 ; % -> the noise power at the input of the PA in dBm
n_rx_dBm = 100 ; % -> the noise power at Rx in dBm
NF = 5 ; % -> noise figure of the amplifier
fc = 28 ; % -> carrier frequency in GHz
K = 5 ; % -> Rician factor
d = 50 ; % -> the distance between Tx and Rx in meters
d_v = 5 ; % -> the vertical distance between Tx and RIS in meters
d_h = 5 ; % -> the horizontal distance between Tx and RIS in meters

[scale_param, shape_param] = find_gamma_parameters(L, N, P_max_dBm, ...
    P_t_dBm, G_max_dB, n_amp_dBm, n_rx_dBm, NF, fc, K, d, d_v, d_h) ;

function [scale_param, shape_param] = find_gamma_parameters(L, N, P_max_dBm, P_t_dBm, G_max_dB, n_amp_dBm, n_rx_dBm, NF, fc, K, d, d_v, d_h)
%{

Copyright (c) 2020-2021 Koc University - CoreLab. 
We kindly ask you to cite our paper named "An Amplifying RIS Architecture 
with a Single Power Amplifier: Energy Efficiency and Error Performance 
Analyses" if you use this code in your studies.


INPUTS:
L -> the number of SNR samples
N -> the RIS size
P_max_dBm -> the maximum power of the PA at the RIS in dBm
P_t_dBm -> transmit power at the BS in dBm
G_max_dB -> maximum gain of the PA in dB
n_amp_dBm -> the noise power at the input of the PA in dBm
n_rx_dBm -> the noise power at Rx in dBm
NF -> noise figure of the amplifier
fc -> carrier frequency in GHz
K -> Rician factor
d -> the distance between Tx and Rx in meters
d_v -> the vertical distance between Tx and RIS in meters
d_h -> the horizontal distance between Tx and RIS in meters

OUTPUTS:
scale_param -> the scale parameter of the estimated Gamma distribution
shape_param -> the shape parameter of the estimated Gamma distribution


%}

d_1 = sqrt(d_h^2 + d_v^2) ; % the distance between Tx and RIS in meters
d_RIS_Rx_hrz = d - d_h ; % the horizontal distance between RIS and Rx in meters
d_2 = sqrt(d_RIS_Rx_hrz^2 + d_v^2) ; % the distance between RIS and Rx in meters

F = 10^(NF / 10) ; % noise factor

sigma_rx = sqrt(10^(n_rx_dBm / 10) * 1e-3) ; % dBm to W
sigma_amp = sqrt(10^(n_amp_dBm / 10) * 1e-3) ; % dBm to W
P_t = 10^(P_t_dBm / 10) * 1e-3 ; % dBm to W

% LOS communication is considered for h
% where NLOS communication is considered for g
lambda_1_dB_LOS = 32.4 + 17.3 * log10(d_1) + 20 * log10(fc) ; % path-loss for h (LOS)
lambda_2_dB_LOS = 32.4 + 17.3 * log10(d_2) + 20 * log10(fc) ; % path-loss for g (LOS)
lambda_2_dB_NLOS = max([32.4 + 20 * log10(fc) + 31.9 * log10(d_2) lambda_2_dB_LOS]) ; % path-loss for g (NLOS)

h_1 = (1 ./ sqrt(10.^(lambda_1_dB_LOS/10))) .* (sqrt(K/(K+1)) + sqrt(1/(K+1)) * 1 / sqrt(2) * ( randn(L, N) + 1i*randn(L, N) ) ) ; % generating the channel h
h_2 = (1 ./ sqrt(10.^(lambda_2_dB_NLOS/10))) .* 1 ./ sqrt(2) .* ( randn(L, N) + 1i*randn(L, N) ) ; % generating the channel g

phi = exp(-1i * angle(h_1) ) ; % phase shift vector of RIS_1
theta = exp(-1i * angle(h_2) ) ; % phase shift vector of RIS_2

G_opt_bar_dB = P_max_dBm - 10*log10(((sum(abs(h_1),2)).^2) * P_t / 0.001) ; % G_opt_bar in dB
G_opt = 10.^( min( ones(L,1) * G_max_dB , G_opt_bar_dB ) / 10) ; % optimized G (G_opt)

% defining the maximized received SNR
h_amp = sqrt(P_t) * sum(sum(h_1 .* phi,2) .* sqrt(G_opt) .* sqrt(1 / N) .* theta .* h_2,2) ;

gamma_act = ( abs(h_amp).^2 ) ...
    ./ (sigma_rx^2 + abs( sqrt(G_opt) / sqrt(N) * sqrt(F) * sigma_amp .* sum(theta .* h_2, 2) ).^2 ) ;

dist_gamma = fitdist(gamma_act,'Gamma') ; % fitting the gamma distribution to the SNR samples
scale_param = dist_gamma.b ; % scale parameter of the fitted gamma distribution
shape_param = dist_gamma.a ; % shape parameter of the fitted gamma distribution
end
