%% verify_presentation_metrics.m
% Thin runner for canonical RMSE and Pd/BER experiment scaffolds.
clear; close all; clc;
addpath('config');
addpath('core/tx'); addpath('core/channel'); addpath('core/rx');
addpath('core/equalization'); addpath('core/sensing'); addpath('core/metrics');
addpath('experiments');

cfg = default_config();
rmse_results = run_rmse_vs_snr(cfg);
pd_ber_results = run_pd_vs_ber(cfg);
save('presentation_metrics_results.mat', 'rmse_results', 'pd_ber_results');
disp('Presentation metrics experiments complete.');
