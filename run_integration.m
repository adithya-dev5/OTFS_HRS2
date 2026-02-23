% run_integration.m
% Thin orchestration runner using canonical config/core pipeline.

clear; close all; clc;
addpath('config');
addpath('core/tx'); addpath('core/channel'); addpath('core/rx');
addpath('core/equalization'); addpath('core/sensing'); addpath('core/metrics'); addpath('core/utils');

cfg = default_config();
[cfg, manifest] = init_experiment_context(cfg);

[dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);

DataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
tx_pkg = build_tx_waveform(cfg, DataSyms);
X = tx_pkg.X;
tx_time = tx_pkg.tx_time;
rx = channel_model(tx_time, cfg, cfg.channel.mode);

rx_power = mean(abs(rx).^2);
noise_var = rx_power / (10^(cfg.channel.noise_snr_db/10));
rx = rx + sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));

assert(exist('tx_time','var')==1 && ~isempty(tx_time), 'Missing canonical tx_time before demodulation.');
Y = otfs_demodulate(rx, cfg);
Heff = build_heff(cfg);
Xhat = mmse_equalize(Y, Heff, noise_var, cfg);
RVM = compute_rvm(Y, cfg, 4, 4);
peaks = find_rvm_peaks(RVM, numel(cfg.targets), cfg);
metrics = evaluate_isac(X, Xhat, peaks, cfg);

save('receiver_results.mat', 'X', 'tx_time', 'Y', 'Xhat', 'RVM', 'peaks', 'metrics', 'cfg', 'manifest');
fprintf('Integration run complete. Results saved to receiver_results.mat\n');
