function [cfg, manifest] = init_experiment_context(cfg)
% INIT_EXPERIMENT_CONTEXT Deterministic init and manifest generation.
assert(isfield(cfg,'seed'), 'cfg.seed is required.');
rng(cfg.seed);
cfg.seed_initialized = true;

if ~isfield(cfg,'git_commit_hash') || isempty(cfg.git_commit_hash)
    [st, out] = system('git rev-parse --short HEAD');
    if st == 0
        cfg.git_commit_hash = strtrim(out);
    else
        cfg.git_commit_hash = 'unknown';
    end
end

validate_physical_units(cfg);

manifest = struct();
manifest.seed = cfg.seed;
manifest.timestamp = datestr(now, 30);
manifest.git_commit_hash = cfg.git_commit_hash;
manifest.cfg = cfg;
end
