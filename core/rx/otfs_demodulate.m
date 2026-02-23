function Y = otfs_demodulate(rx_cp, cfg)
% OTFS_DEMODULATE Preserve existing demodulation mapping behavior.
assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
N = cfg.N; M = cfg.M; Mcp = cfg.Mcp;

rx_cp = rx_cp(:);
if Mcp > 0
    rx_noCP = rx_cp(Mcp+1:end);
else
    rx_noCP = rx_cp;
end

len = N * M;
if numel(rx_noCP) < len
    rx_noCP = [rx_noCP; zeros(len - numel(rx_noCP),1)];
elseif numel(rx_noCP) > len
    rx_noCP = rx_noCP(1:len);
end

X_fast_slow = reshape(rx_noCP, [M, N]);
X_time_delay = X_fast_slow.';
Y = fft(X_time_delay, [], 1);
end
