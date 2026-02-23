function test_channel_scaling()
addpath('config'); addpath('core/channel');
cfg = default_config();
cfg.channel.mode = "synthetic_fractional";
cfg.targets = struct('tau', 0, 'nu', 0, 'gain', 2.0);

tx = ones(cfg.M*cfg.N,1);
rx = channel_model(tx, cfg, cfg.channel.mode);

ratio = mean(abs(rx)) / mean(abs(tx));
assert(abs(ratio - 2.0) < 1e-6, 'Channel gain scaling mismatch.');
end
