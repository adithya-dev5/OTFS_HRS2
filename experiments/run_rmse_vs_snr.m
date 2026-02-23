function results = run_rmse_vs_snr(cfg)
% RUN_RMSE_VS_SNR Canonical RMSE sweep runner (no inline plotting).
if nargin < 1 || isempty(cfg), cfg = default_config(); end
rng(cfg.seed);
if ~exist(cfg.experiments.output_dir, 'dir'), mkdir(cfg.experiments.output_dir); end

Heff = build_heff(cfg);
SNR_dB_Range = cfg.experiments.snr_db_range;
rmse_r = nan(size(SNR_dB_Range));
rmse_v = nan(size(SNR_dB_Range));

for i = 1:length(SNR_dB_Range)
    snr = SNR_dB_Range(i);
    dataSyms = qammod(randi([0 3], cfg.M*cfg.N - 1, 1), 4, 'UnitAveragePower', true);
    [dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
    maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
    maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
    maskGuard = maskGuard & ~maskPilot;
    maskData = ~(maskGuard | maskPilot);
    dataSyms = dataSyms(1:nnz(maskData));

    X = build_dd_grid(cfg, dataSyms);
    tx = otfs_modulate(X, cfg);
    rx_clean = channel_model(tx, cfg, cfg.channel.mode);
    sig_pwr = mean(abs(rx_clean).^2);
    noise_pwr = sig_pwr / (10^(snr/10));
    y_noisy = rx_clean + sqrt(noise_pwr/2)*(randn(size(rx_clean))+1j*randn(size(rx_clean)));

    Y = otfs_demodulate(y_noisy, cfg);
    Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg);
    RVM = compute_rvm(Y, cfg, 4, 4);
    peaks = find_rvm_peaks(RVM, numel(cfg.targets), cfg);
    m = evaluate_isac(X, Xhat, peaks, cfg);

    rmse_r(i) = m.RMSE_range_idx;
    rmse_v(i) = m.RMSE_doppler_idx;
end

results = struct('cfg', cfg, 'snr_db', SNR_dB_Range, 'rmse_range_idx', rmse_r, 'rmse_doppler_idx', rmse_v, 'timestamp', datestr(now,30));
save(fullfile(cfg.experiments.output_dir, 'rmse_vs_snr.mat'), 'results');
end
