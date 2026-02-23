function X = build_dd_grid(cfg, DataSymbols)
% BUILD_DD_GRID Build NxM DD grid using canonical cfg.

assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
M = cfg.M; N = cfg.N; lp = cfg.lp; kp = cfg.kp; lTau = cfg.lTau; kNu = cfg.kNu; Xp = cfg.Xp;

if lp < 0 || lp >= M, error('lp must be 0-based in [0,%d].', M-1); end
if kp < 0 || kp >= N, error('kp must be 0-based in [0,%d].', N-1); end

X = zeros(N, M);
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskPilot = (dd_L == lp) & (dd_K == kp);
maskGuard = (abs(dd_L - lp) <= lTau) | (abs(dd_K - kp) <= 2*kNu);
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);

if numel(DataSymbols) ~= nnz(maskData)
    error('DataSymbols length mismatch: got %d expected %d.', numel(DataSymbols), nnz(maskData));
end

X(kp+1, lp+1) = Xp;
X(maskData) = DataSymbols(:);
end
