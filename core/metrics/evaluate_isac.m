function metrics = evaluate_isac(X, Xhat, peaks, cfg)
% EVALUATE_ISAC Canonical metrics aggregation with no duplicated blocks.
metrics = struct('SER', NaN, 'peaks', peaks, 'RMSE_range_idx', NaN, 'RMSE_doppler_idx', NaN);

if ~isempty(X)
    [dd_L, dd_K] = meshgrid(0:cfg.M-1, 0:cfg.N-1);
    maskPilot = (dd_L == cfg.lp) & (dd_K == cfg.kp);
    maskGuard = (abs(dd_L - cfg.lp) <= cfg.lTau) | (abs(dd_K - cfg.kp) <= 2*cfg.kNu);
    maskGuard = maskGuard & ~maskPilot;
    maskData = ~(maskGuard | maskPilot);

    txsym = X(maskData);
    rxsym = Xhat(maskData);
    tx_dec = sign(real(txsym)) + 1j*sign(imag(txsym));
    rx_dec = sign(real(rxsym)) + 1j*sign(imag(rxsym));
    metrics.SER = sum(tx_dec ~= rx_dec) / max(numel(tx_dec),1);
end

rmse = rmse_from_peaks(peaks, cfg.targets, cfg);
metrics.RMSE_range_idx = rmse.RMSE_range_idx;
metrics.RMSE_doppler_idx = rmse.RMSE_doppler_idx;
end
