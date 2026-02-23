function results = run_pd_vs_ber(cfg)
% RUN_PD_VS_BER Canonical Pd vs BER scaffold runner (no inline plotting).
if nargin < 1 || isempty(cfg), cfg = default_config(); end
[cfg, manifest] = init_experiment_context(cfg);
invariance_report = validate_numerical_invariance(cfg);
if ~exist(cfg.experiments.output_dir, 'dir'), mkdir(cfg.experiments.output_dir); end

if cfg.enable_profiling
    t_heff = tic; Heff = build_heff(cfg); heff_s = toc(t_heff);
else
    Heff = build_heff(cfg); heff_s = NaN;
end
SNR_dB_Range = cfg.experiments.snr_db_range;
ber = nan(size(SNR_dB_Range)); pd = nan(size(SNR_dB_Range));
profile = struct('heff_s', heff_s, 'channel_s', zeros(size(SNR_dB_Range)), ...
    'equalize_s', zeros(size(SNR_dB_Range)), 'rvm_s', zeros(size(SNR_dB_Range)));

for i = 1:length(SNR_dB_Range)
    snr = SNR_dB_Range(i);
    [dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
    maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
    maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
    maskGuard = maskGuard & ~maskPilot;
    maskData = ~(maskGuard | maskPilot);

    dataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
    tx_pkg = build_tx_waveform(cfg, dataSyms);
    X = tx_pkg.X; tx_time = tx_pkg.tx_time;

    if cfg.enable_profiling
        tic; rx_clean = channel_model(tx_time, cfg, cfg.channel.mode); profile.channel_s(i)=toc;
    else
        rx_clean = channel_model(tx_time, cfg, cfg.channel.mode);
    end

    sig_pwr = mean(abs(rx_clean).^2);
    noise_pwr = sig_pwr / (10^(snr/10));
    y_noisy = rx_clean + sqrt(noise_pwr/2)*(randn(size(rx_clean))+1j*randn(size(rx_clean)));

    assert(exist('tx_time','var')==1 && ~isempty(tx_time), 'Missing canonical tx_time before demodulation.');
    Y = otfs_demodulate(y_noisy, cfg);

    if cfg.enable_profiling
        tic; Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg); profile.equalize_s(i)=toc;
    else
        Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg);
    end
    [ber(i), ~, ~] = ber_qpsk(X(maskData), Xhat(maskData));

    if cfg.enable_profiling
        tic; RVM = compute_rvm(Y, cfg, 4, 4); profile.rvm_s(i)=toc;
    else
        RVM = compute_rvm(Y, cfg, 4, 4);
    end
    peaks = find_rvm_peaks(RVM, numel(cfg.targets), cfg);
    pd(i) = double(~isempty(peaks));
end

results = struct('cfg', cfg, 'seed', cfg.seed, 'snr_db', SNR_dB_Range, ...
    'pd', pd, 'ber', ber, 'timestamp', manifest.timestamp, ...
    'git_commit_hash', manifest.git_commit_hash, 'profile', profile, 'invariance_report', invariance_report);
save(fullfile(cfg.experiments.output_dir, 'pd_vs_ber.mat'), 'results', 'manifest');
end
