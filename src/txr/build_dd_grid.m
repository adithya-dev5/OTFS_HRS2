function X = build_dd_grid(M, N, lp, kp, lTau, kNu, DataSymbols, Xp)
%% BUILD_DD_GRID Construct N×M delay-Doppler grid with pilot + guards + data.
%  X = build_dd_grid(M,N,lp,kp,lTau,kNu,DataSymbols,Xp)

% convert to doubles
M=double(M); N=double(N); lp=double(lp); kp=double(kp); lTau=double(lTau); kNu=double(kNu);
% heuristically accept 1-based indices: if lp in 1..M assume 1-based
% --- FIX: Strict 0-based Indexing ---
if lp < 0 || lp >= M, error('Pilot Delay Index (lp) must be 0-based (0 to %d). Got: %d', M-1, lp); end
if kp < 0 || kp >= N, error('Pilot Doppler Index (kp) must be 0-based (0 to %d). Got: %d', N-1, kp); end% initialize
X = zeros(N, M);
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskPilot = (dd_L == lp) & (dd_K == kp);
maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);
numData = nnz(maskData);
if numel(DataSymbols) ~= numData
    error('DataSymbols length (%d) does not match available data slots (%d).', numel(DataSymbols), numData);
end
% place pilot
X(kp+1, lp+1) = Xp;
% place data in column-major order consistent with MATLAB linear indexing
dataIdx = find(maskData);
X(dataIdx) = DataSymbols(:);
end
