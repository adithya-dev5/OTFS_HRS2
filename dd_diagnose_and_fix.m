% dd_diagnose_and_fix.m
% Run from project root. Requires waveform_package.mat and receiver_results.mat

clear; close all; clc;
wp = load('waveform_package.mat');
rr = load('receiver_results.mat');

% get params, X (tx DD) and candidate rx domains
params = [];
if isfield(wp,'params'), params = wp.params; end
if isfield(wp,'X'), Xtx = wp.X; else Xtx = []; end

% get candidate DD from receiver_results
Y_rx_candidates = struct();

if isfield(rr,'Y')
    Y_rx_candidates.Y_from_rr = rr.Y;
end
if isfield(rr,'Xhat')
    Y_rx_candidates.Xhat = rr.Xhat;   % possibly post-eq
end

% get time-domain tx and rx if present
tx = [];
rx = [];
if isfield(wp,'tx_signal'), tx = wp.tx_signal(:); elseif isfield(wp,'tx'), tx = wp.tx(:); end
if isfield(rr,'rx'), rx = rr.rx(:); elseif isfield(rr,'rx_signal'), rx = rr.rx_signal(:); end

% attempt to infer M,N
if isfield(params,'M') && isfield(params,'N')
    M = params.M; N = params.N;
else
    % attempt to infer from Xtx
    if ~isempty(Xtx)
        [r,c] = size(Xtx);
        % guess which dimension is which: prefer N x M
        if r <= c
            N = r; M = c;
        else
            M = r; N = c;
        end
    else
        error('Cannot infer M,N. Provide params.M and params.N or X in waveform_package.mat');
    end
end
fprintf('Using M=%d, N=%d\n',M,N);

% helper demod function (two common conventions)
demod_try = @(sig, Mcp, Mloc, Nloc, mode) ...
    demod_helper(sig, Mcp, Mloc, Nloc, mode);

% 1) If we have time-domain rx, produce two demod options
modes = {'fft_dim1_then_dim2','fft_dim2_then_dim1','ifftshift_option1','ifftshift_option2'};
Y_from_rx_variants = struct();
if ~isempty(rx)
    Mcp = 0; if isfield(params,'Mcp'), Mcp = params.Mcp; end
    for im=1:numel(modes)
        mode = modes{im};
        Yv = demod_try(rx, Mcp, M, N, mode);
        Y_from_rx_variants.(mode) = Yv;
    end
    Y_rx_candidates.from_rx = Y_from_rx_variants;
end

% 2) For each candidate Y, compute pilot energy and compare with X
% Determine pilot positions from Xtx (if available)
if ~isempty(Xtx)
    % try to map Xtx into N x M convention used by code
    [rX,cX] = size(Xtx);
    if rX==N && cX==M
        Xtx_dd = Xtx;
    elseif rX==M && cX==N
        Xtx_dd = Xtx.'; 
    else
        % ambiguous: reshape
        Xtx_dd = reshape(Xtx, [N, M]);
    end
    % pilot mask: non-zero entries in Xtx
    maskPilot = (abs(Xtx_dd) > 0);
    pilot_coords = find(maskPilot);
    [pil_r, pil_c] = ind2sub(size(Xtx_dd), pilot_coords);
    fprintf('TX pilot/guard nonzero count: %d\n', numel(pilot_coords));
else
    Xtx_dd = [];
    pilot_coords = [];
    pil_r=[]; pil_c=[];
    fprintf('No X in waveform_package.mat. Skipping pilot-location tests.\n');
end

% function to score a candidate Y: pilot energy inside pilot mask vs total energy
score_candidate = @(Y) score_Y_against_pilot(Y, Xtx_dd);

fprintf('\nCandidate checks and scores:\n');
candidate_names = fieldnames(Y_rx_candidates);
scores = struct();
for k=1:numel(candidate_names)
    key = candidate_names{k};
    val = Y_rx_candidates.(key);
    if isstruct(val)  % variants
        subnames = fieldnames(val);
        for s=1:numel(subnames)
            name_full = sprintf('%s.%s',key,subnames{s});
            Ycand = val.(subnames{s});
            sc = score_candidate(Ycand);
            scores.(name_full) = sc;
            fprintf('%-30s : pilot_energy_fraction = %.4f   (max mags: tx=%.3g, rx=%.3g)\n', name_full, sc.frac_pilot, sc.max_tx, sc.max_rx);
        end
    else
        Ycand = val;
        sc = score_candidate(Ycand);
        scores.(key) = sc;
        fprintf('%-30s : pilot_energy_fraction = %.4f   (max mags: tx=%.3g, rx=%.3g)\n', key, sc.frac_pilot, sc.max_tx, sc.max_rx);
    end
end

% Show the best candidate by frac_pilot
allnames = fieldnames(scores);
best = allnames{1}; bestval = scores.(best).frac_pilot;
for i=2:numel(allnames)
    if scores.(allnames{i}).frac_pilot > bestval
        best = allnames{i}; bestval = scores.(allnames{i}).frac_pilot;
    end
end
fprintf('\nBest candidate mapping appears to be: %s (pilot fraction = %.4f)\n', best, bestval);

% extract best Y and plot TX vs best RX DD
if contains(best, '.')
    tmp = strsplit(best,'.');
    base = tmp{1}; sub = tmp{2};
    Ybest = Y_rx_candidates.(base).(sub);
else
    Ybest = Y_rx_candidates.(best);
end

% Normalize and plot
figure; imagesc(abs(Xtx_dd)/max(abs(Xtx_dd(:))+eps)); axis xy; colorbar; title('TX (normalized)'); xlabel('Delay'); ylabel('Doppler');
figure; imagesc(abs(Ybest)/max(abs(Ybest(:))+eps)); axis xy; colorbar; title(['Best RX candidate: ' best ' (normalized)']); xlabel('Delay'); ylabel('Doppler');

% If pilot exists in time-domain tx, cross-correlate with rx to see pilot energy
if ~isempty(tx) && ~isempty(rx)
    % find the pilot waveform in time-domain: if tx is long, we search for high-energy segment
    [~, idx_max] = max(abs(tx));
    pilot_time_sample = idx_max;
    % compute cross-correlation of rx with pilot time waveform (use full tx as probe)
    xc = xcorr(rx, tx, 'coeff');
    [~, ix] = max(abs(xc));
    lag = (ix - (numel(tx)));
    fprintf('Cross-correlation peak between rx and tx at lag = %d samples (positive means rx delayed)\n', lag);
    figure; plot(xc); title('xcorr(rx,tx)'); grid on;
end

% Helper functions -------------------------------------------------------
function Y = demod_helper(sig, Mcp, Mloc, Nloc, mode)
    s = sig(:);
    if Mcp > 0 && numel(s) > Mcp, s = s(Mcp+1:end); end
    len = Mloc * Nloc;
    if numel(s) < len, s = [s; zeros(len-numel(s),1)]; end
    s = s(1:len);
    switch mode
        case 'fft_dim1_then_dim2'
            tf = reshape(s, [Nloc, Mloc]);      % N x M
            Y = fft(fft(tf, [], 1), [], 2);
        case 'fft_dim2_then_dim1'
            tf = reshape(s, [Mloc, Nloc]);      % M x N
            Ytmp = fft(fft(tf, [], 2), [], 1);
            Y = Ytmp.'; % bring back to N x M
        case 'ifftshift_option1'
            tf = reshape(s, [Nloc, Mloc]);
            tf = ifftshift(tf,1);
            Y = fft(fft(tf, [], 1), [], 2);
        case 'ifftshift_option2'
            tf = reshape(s, [Nloc, Mloc]);
            tf = ifftshift(tf,2);
            Y = fft(fft(tf, [], 1), [], 2);
        case 'true_otfs'
        % The standard OTFS: FFT along Doppler (1), IFFT along Delay (2)
        tf = reshape(s, [Nloc, Mloc]);
        Y = ifft(fft(tf, [], 1), [], 2);
        otherwise
            tf = reshape(s, [Nloc, Mloc]);
            Y = fft(fft(tf, [], 1), [], 2);
    end
end

function sc = score_Y_against_pilot(Y, Xtx_dd)
    sc = struct('frac_pilot',0,'sum_pilot',0,'sum_total',0,'max_tx',0,'max_rx',0);
    magY = abs(Y);
    sc.max_rx = max(magY(:));
    if ~isempty(Xtx_dd)
        magX = abs(Xtx_dd);
        sc.max_tx = max(magX(:));
        mask = (magX > (0.5*max(magX(:)))); % strong pilot area (threshold)
        if ~any(mask(:))
            % fallback to any nonzero
            mask = (magX > 0);
        end
        sc.sum_pilot = sum(magY(mask));
        sc.sum_total = sum(magY(:)) + eps;
        sc.frac_pilot = sc.sum_pilot / sc.sum_total;
    else
        sc.frac_pilot = 0;
        sc.sum_total = sum(magY(:));
        sc.sum_pilot = 0;
        sc.max_tx = NaN;
    end
end
