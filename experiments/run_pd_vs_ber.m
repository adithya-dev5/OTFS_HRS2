function results = run_pd_vs_ber(cfg)
% RUN_PD_VS_BER Canonical Pd vs BER scaffold runner (no inline plotting).
if nargin < 1 || isempty(cfg), cfg = default_config(); end
rng(cfg.seed);
if ~exist(cfg.experiments.output_dir, 'dir'), mkdir(cfg.experiments.output_dir); end

Heff = build_heff(cfg);
SNR_dB_Range = cfg.experiments.snr_db_range;
ber = nan(size(SNR_dB_Range));
pd = nan(size(SNR_dB_Range));

for i = 1:length(SNR_dB_Range)
    snr = SNR_dB_Range(i);
    [dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
    maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
    maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
    maskGuard = maskGuard & ~maskPilot;
    maskData = ~(maskGuard | maskPilot);

    dataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
    X = build_dd_grid(cfg, dataSyms);
    tx = otfs_modulate(X, cfg);
    rx_clean = channel_model(tx, cfg, cfg.channel.mode);

    sig_pwr = mean(abs(rx_clean).^2);
    noise_pwr = sig_pwr / (10^(snr/10));
    Y = otfs_demodulate(rx_clean + sqrt(noise_pwr/2)*(randn(size(rx_clean))+1j*randn(size(rx_clean))), cfg);

    Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg);
    [ber(i), ~, ~] = ber_qpsk(X(maskData), Xhat(maskData));

    RVM = compute_rvm(Y, cfg, 4, 4);
    peaks = find_rvm_peaks(RVM, numel(cfg.targets), cfg);
    pd(i) = double(~isempty(peaks)); % structural placeholder (binary detect-any)
end

results = struct('cfg', cfg, 'snr_db', SNR_dB_Range, 'pd', pd, 'ber', ber, 'timestamp', datestr(now,30));
save(fullfile(cfg.experiments.output_dir, 'pd_vs_ber.mat'), 'results');
end
