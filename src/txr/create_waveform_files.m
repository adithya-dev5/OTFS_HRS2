% create_waveform_files.m
% Run this once to generate the helper .m files for the OTFS waveform + channel utilities.

files = {};

files{end+1} = struct('name','build_dd_grid.m','content',[
"function X = build_dd_grid(M, N, lp, kp, lTau, kNu, DataSymbols, Xp)"
"%% BUILD_DD_GRID Construct N×M delay-Doppler grid with pilot + guards + data."
"%  X = build_dd_grid(M,N,lp,kp,lTau,kNu,DataSymbols,Xp)"
""
"% convert to doubles"
"M=double(M); N=double(N); lp=double(lp); kp=double(kp); lTau=double(lTau); kNu=double(kNu);"
"% heuristically accept 1-based indices: if lp in 1..M assume 1-based"
"if lp>=1 && lp<=M && lp==floor(lp) && lp~=0; lp = lp-1; end"
"if kp>=1 && kp<=N && kp==floor(kp) && kp~=0; kp = kp-1; end"
"% initialize"
"X = zeros(N, M);"
"[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);"
"maskPilot = (dd_L == lp) & (dd_K == kp);"
"maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );"
"maskGuard = maskGuard & ~maskPilot;"
"maskData = ~(maskGuard | maskPilot);"
"numData = nnz(maskData);"
"if numel(DataSymbols) ~= numData"
"    error('DataSymbols length (%d) does not match available data slots (%d).', numel(DataSymbols), numData);"
"end"
"% place pilot"
"X(kp+1, lp+1) = Xp;"
"% place data in column-major order consistent with MATLAB linear indexing"
"dataIdx = find(maskData);"
"X(dataIdx) = DataSymbols(:);"
"end"
]);

files{end+1} = struct('name','otfs_modulate.m','content',[
"function tx = otfs_modulate(X, N, M)"
"% OTFS modulate: IFFT across rows then reshape to time-domain vector."
"[Na, Ma] = size(X);"
"assert(Na==N && Ma==M,'X size mismatch with (N,M)');"
"% IFFT across rows (dimension 1)"
"timeCols = ifft(X, [], 1);"
"% reshape columns into N*M vector (column-major)"
"tx = reshape(timeCols, [], 1);"
"% normalize to unit average power"
"p = mean(abs(tx).^2);"
"if p>0; tx = tx / sqrt(p); end"
"end"
]);

files{end+1} = struct('name','apply_channel.m','content',[
"function rx = apply_channel(tx, targets, params)"
"% APPLY_CHANNEL: fractional-delay (freq-domain) + time-domain Doppler per path."
"tx = tx(:);"
"if ~isfield(params,'fs'), error('params.fs is required'); end"
"fs = params.fs;"
"if ~isfield(params,'Lfft') || isempty(params.Lfft)"
"    Lfft = 2^nextpow2(max(4*length(tx), length(tx)));"
"else"
"    Lfft = max(params.Lfft, length(tx));"
"end"
"t0 = 0; if isfield(params,'t0'), t0 = params.t0; end"
"tvec = t0 + (0:length(tx)-1).' / fs;"
"fvec = (0:Lfft-1).' * (fs / Lfft);"
"TX_F = fft(tx, Lfft);"
"rx = zeros(size(tx));"
"for p = 1:numel(targets)"
"    tg = targets(p);"
"    if ~isfield(tg,'tau') || ~isfield(tg,'nu') || ~isfield(tg,'gain'), error('target missing fields'); end"
"    tau = tg.tau; nu = tg.nu; gain = tg.gain;"
"    phase_delay = exp(-1j*2*pi * fvec * tau);"
"    X_shifted = TX_F .* phase_delay;"
"    x_delayed_full = ifft(X_shifted, Lfft);"
"    x_delayed = x_delayed_full(1:length(tx));"
"    doppler_phasor = exp(1j*2*pi * nu * tvec);"
"    x_doppler = x_delayed .* doppler_phasor;"
"    rx = rx + gain * x_doppler;"
"end"
"end"
]);

files{end+1} = struct('name','add_cp.m','content',[
"function xcp = add_cp(x, Mcp)"
"x = x(:);"
"if Mcp==0, xcp = x; return; end"
"L = length(x); assert(Mcp < L, 'Mcp must be less than signal length');"
"prefix = x(end-Mcp+1:end);"
"xcp = [prefix; x];"
"end"
]);

files{end+1} = struct('name','save_tx_for_receiver.m','content',[
"function save_tx_for_receiver(filename, tx_struct)"
"if ~(ischar(filename) || isstring(filename)), error('filename must be string'); end"
"if ~isfield(tx_struct,'tx_signal') || ~isfield(tx_struct,'X') || ~isfield(tx_struct,'params')"
"    error('tx_struct must contain tx_signal, X, and params');"
"end"
"tx_struct.tx_signal = tx_struct.tx_signal(:);"
"save(filename, '-struct', 'tx_struct');"
"fprintf('Saved %s (tx len=%d)\\n', filename, length(tx_struct.tx_signal));"
"end"
]);

files{end+1} = struct('name','waveform_and_channel.m','content',[
"function [X, tx, params, targets] = waveform_and_channel(cfg)"
"if ~exist('cfg','var') || isempty(cfg), cfg = struct(); end"
"defaults = struct('M',32,'N',32,'lp',floor(32/2)-1,'kp',floor(32/2)-1,'lTau',2,'kNu',2,'Xp',1+0j,'fs',15.36e6,'Mcp',64);"
"flds = fieldnames(defaults);"
"for i=1:numel(flds), if ~isfield(cfg,flds{i}) || isempty(cfg.(flds{i})), cfg.(flds{i}) = defaults.(flds{i}); end, end"
"M = cfg.M; N = cfg.N; lp = cfg.lp; kp = cfg.kp; lTau = cfg.lTau; kNu = cfg.kNu; Xp = cfg.Xp;"
"[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);"
"maskPilot = (dd_L == lp) & (dd_K == kp);"
"maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );"
"maskGuard = maskGuard & ~maskPilot;"
"maskData = ~(maskGuard | maskPilot);"
"numData = nnz(maskData);"
"DataSymbols = (1-2*randi([0 1], numData,1)) + 1j*(1-2*randi([0 1], numData,1));"
"DataSymbols = DataSymbols / sqrt(mean(abs(DataSymbols).^2));"
"X = build_dd_grid(M, N, lp, kp, lTau, kNu, DataSymbols, Xp);"
"tx = otfs_modulate(X, N, M);"
"params.fs = cfg.fs; params.M = M; params.N = N; params.Mcp = cfg.Mcp;"
"params.Lfft = 2^nextpow2(4*length(tx));"
"params.fvec = (0:params.Lfft-1).' * (params.fs / params.Lfft);"
"params.tvec = (0:length(tx)-1).' / params.fs;"
"if isfield(cfg,'targets'), targets = cfg.targets; else targets = []; end"
"end"
]);

files{end+1} = struct('name','package_and_save.m','content',[
"%% package_and_save.m"
"cfg = struct('M',32,'N',32,'lp',floor(32/2)-1,'kp',floor(32/2)-1,'lTau',2,'kNu',2,'Xp',1+0j,'fs',15.36e6,'Mcp',64);"
"[X, tx, params, ~] = waveform_and_channel(cfg);"
"targets(1).tau = 1.2e-6; targets(1).nu = 3; targets(1).gain = 0.8*exp(1j*0.2);"
"targets(2).tau = 3.8e-6; targets(2).nu = -1.5; targets(2).gain = 0.4*exp(-1j*0.5);"
"waveform_package.X = X; waveform_package.tx_signal = tx; waveform_package.params = params;"
"waveform_package.targets = targets; waveform_package.Lfft = params.Lfft;"
"waveform_package.fvec = params.fvec; waveform_package.tvec = params.tvec;"
"save('waveform_package.mat','-struct','waveform_package');"
"fprintf('Saved waveform_package.mat into current folder.\\n');"
]);

% write files
for k=1:numel(files)
    fid = fopen(files{k}.name,'w');
    if fid==-1, error('Cannot open file %s for writing', files{k}.name); end
    fprintf(fid, '%s\n', files{k}.content);
    fclose(fid);
    fprintf('Wrote %s\n', files{k}.name);
end

fprintf('\\nDone. Run: addpath(pwd); then try: [X,tx,params,targets]=waveform_and_channel();\\n');
