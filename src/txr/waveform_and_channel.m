function [X, tx, params, targets] = waveform_and_channel(cfg)
if ~exist('cfg','var') || isempty(cfg), cfg = struct(); end
defaults = struct('M',32,'N',32,'lp',floor(32/2)-1,'kp',floor(32/2)-1,'lTau',2,'kNu',2,'Xp',1+0j,'fs',15.36e6,'Mcp',64);
flds = fieldnames(defaults);
for i=1:numel(flds), if ~isfield(cfg,flds{i}) || isempty(cfg.(flds{i})), cfg.(flds{i}) = defaults.(flds{i}); end, end
M = cfg.M; N = cfg.N; lp = cfg.lp; kp = cfg.kp; lTau = cfg.lTau; kNu = cfg.kNu; Xp = cfg.Xp;
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskPilot = (dd_L == lp) & (dd_K == kp);
maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);
numData = nnz(maskData);
DataSymbols = (1-2*randi([0 1], numData,1)) + 1j*(1-2*randi([0 1], numData,1));
DataSymbols = DataSymbols / sqrt(mean(abs(DataSymbols).^2));
X = build_dd_grid(M, N, lp, kp, lTau, kNu, DataSymbols, Xp);
tx = otfs_modulate(X, N, M);
params.fs = cfg.fs; params.M = M; params.N = N; params.Mcp = cfg.Mcp;
params.Lfft = 2^nextpow2(4*length(tx));
params.fvec = (0:params.Lfft-1).' * (params.fs / params.Lfft);
params.tvec = (0:length(tx)-1).' / params.fs;
if isfield(cfg,'targets'), targets = cfg.targets; else targets = []; end
end
