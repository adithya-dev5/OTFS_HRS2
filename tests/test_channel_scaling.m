function test_channel_scaling()
% TEST_CHANNEL_SCALING Independent of current working directory.
this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
repo_root = fileparts(this_dir);
addpath(fullfile(repo_root, 'config'));
addpath(fullfile(repo_root, 'core', 'channel'));

cfg = default_config();
cfg.seed_initialized = true; % satisfy channel contract for direct unit test
cfg.channel.mode = "synthetic_fractional";
cfg.targets = struct('tau', 0, 'nu', 0, 'gain', 2.0);

tx = ones(cfg.M*cfg.N,1);
rx = channel_model(tx, cfg, cfg.channel.mode);

ratio = mean(abs(rx)) / mean(abs(tx));
assert(abs(ratio - 2.0) < 1e-6, 'Channel gain scaling mismatch.');
end
