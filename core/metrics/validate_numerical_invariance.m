function report = validate_numerical_invariance(cfg)
% VALIDATE_NUMERICAL_INVARIANCE Compare pre-refactor and post-refactor numerics.
% Halts (error) if deviation exceeds tolerance.

tol = 1e-10;
root = pwd;

pre = run_pre_refactor(cfg, root);
post = run_post_refactor(cfg, root);

report = struct();
report.dd = compare_arrays(pre.X, post.X);
report.channel = compare_arrays(pre.rx, post.rx);
report.rvm = compare_arrays(pre.RVM, post.RVM);
report.ber = compare_scalars(pre.ber, post.ber);
report.rmse = compare_scalars(pre.rmse, post.rmse);
report.tolerance = tol;

fields = {'dd','channel','rvm','ber','rmse'};
for i = 1:numel(fields)
    f = fields{i};
    if report.(f).max_abs_diff > tol || report.(f).max_rel_diff > tol
        error('Numerical invariance failed for %s (abs=%g rel=%g > %g).', ...
            f, report.(f).max_abs_diff, report.(f).max_rel_diff, tol);
    end
end
end

function out = run_pre_refactor(cfg, root)
restoredefaultpath;
addpath(fullfile(root,'src/txr'));
addpath(fullfile(root,'src/receiver'));
rng(cfg.seed);

[dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);
dataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);

X = build_dd_grid(cfg.M, cfg.N, cfg.lp, cfg.kp, cfg.lTau, cfg.kNu, dataSyms, cfg.Xp);
tx = otfs_modulate(X, cfg.N, cfg.M);
params = struct('M',cfg.M,'N',cfg.N,'Mcp',cfg.Mcp,'fs',cfg.fs,'Lfft',cfg.channel.Lfft,'targets',cfg.targets);
rx = apply_channel(tx, cfg.targets, params);

snr = cfg.experiments.snr_db_range(1);
sig_pwr = mean(abs(rx).^2);
noise_pwr = sig_pwr/(10^(snr/10));
noise = sqrt(noise_pwr/2)*(randn(size(rx))+1j*randn(size(rx)));
Y = otfs_demodulate(rx + noise, cfg.Mcp, cfg.N, cfg.M);
Heff = build_Heff(cfg.N, cfg.M, params, []);
Xhat = mmse_equalize(Y, Heff, noise_pwr);
RVM = compute_highres_RVM(Y, 4, 4);

[ber, rmse] = scalar_metrics(X, Xhat, maskData, RVM, cfg);
out = struct('X',X,'rx',rx,'RVM',RVM,'ber',ber,'rmse',rmse);
end

function out = run_post_refactor(cfg, root)
restoredefaultpath;
addpath(fullfile(root,'config'));
addpath(fullfile(root,'core/tx'));
addpath(fullfile(root,'core/channel'));
addpath(fullfile(root,'core/rx'));
addpath(fullfile(root,'core/equalization'));
addpath(fullfile(root,'core/sensing'));
addpath(fullfile(root,'core/metrics'));
addpath(fullfile(root,'core/utils'));

[cfg, ~] = init_experiment_context(cfg);
[dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);
dataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);

pkg = build_tx_waveform(cfg, dataSyms);
X = pkg.X; tx_time = pkg.tx_time;
rx = channel_model(tx_time, cfg, cfg.channel.mode);

snr = cfg.experiments.snr_db_range(1);
sig_pwr = mean(abs(rx).^2);
noise_pwr = sig_pwr/(10^(snr/10));
noise = sqrt(noise_pwr/2)*(randn(size(rx))+1j*randn(size(rx)));
Y = otfs_demodulate(rx + noise, cfg);
Heff = build_heff(cfg);
Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg);
RVM = compute_rvm(Y, cfg, 4, 4);

[ber, rmse] = scalar_metrics(X, Xhat, maskData, RVM, cfg);
out = struct('X',X,'rx',rx,'RVM',RVM,'ber',ber,'rmse',rmse);
end

function [ber, rmse] = scalar_metrics(X, Xhat, maskData, RVM, cfg)
rx_bits = qamdemod(Xhat(maskData), 4, 'UnitAveragePower', true);
tx_bits = qamdemod(X(maskData), 4, 'UnitAveragePower', true);
[~, br] = biterr(tx_bits, rx_bits, 2);
ber = br;

% RMSE over detected target index (first target)
% canonical axis: row=doppler, col=delay
[~, idx] = max(RVM(:));
[r,c] = ind2sub(size(RVM), idx);
q_d = size(RVM,1)/cfg.N;
q_l = size(RVM,2)/cfg.M;
est_ki = r / q_d;
est_li = c / q_l;
true_li = cfg.targets(1).tau * cfg.fs;
true_ki = cfg.targets(1).nu / (cfg.fs/(cfg.M*cfg.N));
rmse = sqrt(mean([(est_li-true_li)^2, (est_ki-true_ki)^2]));
end

function d = compare_arrays(a,b)
absd = abs(a-b);
reld = absd ./ max(abs(a), eps);
d = struct('max_abs_diff', max(absd(:)), 'max_rel_diff', max(reld(:)));
end

function d = compare_scalars(a,b)
absd = abs(a-b);
reld = absd / max(abs(a), eps);
d = struct('max_abs_diff', absd, 'max_rel_diff', reld);
end
