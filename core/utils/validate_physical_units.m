function derived = validate_physical_units(cfg)
% VALIDATE_PHYSICAL_UNITS Enforce physical-unit consistency from cfg.
c = 3e8;

derived = struct();
derived.lambda = c / cfg.fc;
derived.Delta_f = cfg.B / cfg.M;
derived.T = 1 / derived.Delta_f;
derived.frameT = cfg.N * derived.T;

assert(abs(cfg.Delta_f - derived.Delta_f) < 1e-15, 'cfg.Delta_f mismatch.');
assert(abs(cfg.T - derived.T) < 1e-15, 'cfg.T mismatch.');
if isfield(cfg,'Ts')
    assert(abs(cfg.Ts - 1/(cfg.M*derived.Delta_f)) < 1e-15, 'cfg.Ts mismatch.');
end
end
