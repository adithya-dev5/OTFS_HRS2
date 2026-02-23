function metrics = evaluate_ISAC(X, Xhat, params, peaks)
% EVALUATE_ISAC Backward-compatible wrapper to canonical metrics module.
% This wrapper removes duplicated metric logic and delegates to core/metrics/evaluate_isac.

if exist('config/default_config.m', 'file') && exist('core/metrics/evaluate_isac.m', 'file')
    addpath('config');
    addpath('core/metrics');
    cfg = default_config();

    % Preserve caller-supplied parameter overrides where present.
    f = fieldnames(params);
    for i = 1:numel(f)
        cfg.(f{i}) = params.(f{i});
    end

    metrics = evaluate_isac(X, Xhat, peaks, cfg);
else
    % Fallback minimal behavior if core path unavailable.
    metrics = struct('SER', NaN, 'peaks', peaks, 'RMSE_range_idx', NaN, 'RMSE_doppler_idx', NaN);
end
end
