function tx_pkg = build_tx_waveform(cfg, DataSymbols)
% BUILD_TX_WAVEFORM Build DD grid and preserve canonical time-domain waveform.
X = build_dd_grid(cfg, DataSymbols);
tx_time = otfs_modulate(X, cfg);

assert(isvector(tx_time) && ~isempty(tx_time), 'tx_time must exist as non-empty time-domain waveform.');
assert(abs(cfg.fs) > 0, 'cfg.fs must be globally defined and positive.');

tx_pkg = struct('X', X, 'tx_time', tx_time, 'fs', cfg.fs);
end
