function tx = otfs_modulate(X, cfg)
% OTFS_MODULATE Preserve existing modulation mapping behavior.
assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
N = cfg.N; M = cfg.M;

[Na, Ma] = size(X);
if Na ~= N || Ma ~= M
    if Na == M && Ma == N
        X = X.';
    else
        error('X size mismatch. Expected %dx%d.', N, M);
    end
end

X_time_delay = ifft(X, [], 1);
X_fast_slow = X_time_delay.';
tx = X_fast_slow(:);

if isfield(cfg,'normalize') && isfield(cfg.normalize,'tx') && cfg.normalize.tx
    p = mean(abs(tx).^2);
    if p > 0, tx = tx / sqrt(p); end
end
end
