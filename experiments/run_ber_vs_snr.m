function results = run_ber_vs_snr(cfg)
% RUN_BER_VS_SNR Canonical BER sweep runner (no inline plotting).
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
BER = zeros(size(SNR_dB_Range));

profile = struct('heff_s', heff_s, 'channel_s', zeros(size(SNR_dB_Range)), ...
    'equalize_s', zeros(size(SNR_dB_Range)), 'rvm_s', zeros(size(SNR_dB_Range)));

for i = 1:length(SNR_dB_Range)
    snr = SNR_dB_Range(i);
    total_errors = 0; total_bits = 0; frame_count = 0;
    t_ch = 0; t_eq = 0;

    while (total_errors < cfg.experiments.min_errors && frame_count < cfg.experiments.max_frames)
        frame_count = frame_count + 1;

        [dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
        maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
        maskGuard = (abs(dd_L-cfg.lp)<=cfg.lTau) | (abs(dd_K-cfg.kp)<=2*cfg.kNu);
        maskGuard = maskGuard & ~maskPilot;
        maskData = ~(maskGuard | maskPilot);

        dataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
        tx_pkg = build_tx_waveform(cfg, dataSyms);
        X = tx_pkg.X; tx_time = tx_pkg.tx_time;

        if cfg.enable_profiling
            tic; rx_clean = channel_model(tx_time, cfg, cfg.channel.mode); t_ch = t_ch + toc;
        else
            rx_clean = channel_model(tx_time, cfg, cfg.channel.mode);
        end

        sig_pwr = mean(abs(rx_clean).^2);
        noise_pwr = sig_pwr / (10^(snr/10));
        noise = sqrt(noise_pwr/2) * (randn(size(rx_clean)) + 1j*randn(size(rx_clean)));

        assert(exist('tx_time','var')==1 && ~isempty(tx_time), 'Missing canonical tx_time before demodulation.');
        Y = otfs_demodulate(rx_clean + noise, cfg);

        if cfg.enable_profiling
            tic; Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg); t_eq = t_eq + toc;
        else
            Xhat = mmse_equalize(Y, Heff, noise_pwr, cfg);
        end
        if cfg.normalize.rx_pilot_rescale
            pilot_est = Xhat(cfg.kp+1, cfg.lp+1);
            if abs(pilot_est) > 0.01
                Xhat = Xhat * (cfg.Xp / pilot_est);
            end
        end

        [~, errs, bits] = ber_qpsk(X(maskData), Xhat(maskData));
        total_errors = total_errors + errs;
        total_bits = total_bits + bits;
    end

    BER(i) = total_errors / max(total_bits,1);
    profile.channel_s(i) = t_ch;
    profile.equalize_s(i) = t_eq;
end

results = struct();
results.cfg = cfg;
results.seed = cfg.seed;
results.snr_db = SNR_dB_Range;
results.ber = BER;
results.timestamp = manifest.timestamp;
results.git_commit_hash = manifest.git_commit_hash;
results.profile = profile;
results.invariance_report = invariance_report;
save(fullfile(cfg.experiments.output_dir, 'ber_vs_snr.mat'), 'results', 'manifest');
end
