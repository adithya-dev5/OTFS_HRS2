%% generate_ber_curve.m
% Thin runner for canonical BER experiment.
clear; close all; clc;
addpath('config');
addpath('core/tx'); addpath('core/channel'); addpath('core/rx');
addpath('core/equalization'); addpath('core/sensing'); addpath('core/metrics');
addpath('experiments');

cfg = default_config();
results = run_ber_vs_snr(cfg);
save('ber_vs_snr_results.mat', 'results');
disp('BER experiment complete. Saved ber_vs_snr_results.mat');
