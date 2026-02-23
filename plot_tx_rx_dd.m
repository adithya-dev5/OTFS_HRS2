% plot_tx_rx_dd.m
% Visualize transmitted and received signals in Delay-Doppler (DD) domain.
% Requires: waveform_package.mat (from Person A) and receiver_results.mat (from run_integration).
% Saves figures to current folder and prints summary diagnostics.

clear; close all; clc;

% ---- Load data ----
if ~isfile('waveform_package.mat'), error('waveform_package.mat not found'); end
if ~isfile('receiver_results.mat'), error('receiver_results.mat not found'); end

W = load('waveform_package.mat');
R = load('receiver_results.mat');

% Pick M,N and CP from params or X
if isfield(W,'params')
    params = W.params;
else
    params = struct();
end

% decide M/N
if isfield(params,'M') && isfield(params,'N')
    M = params.M; N = params.N;
elseif isfield(W,'X')
    [M,N] = size(W.X); % X is M x N or N x M depending on convention - we try to detect later
else
    error('Cannot determine M and N. Add params.M and params.N to waveform_package.mat or include X.');
end

% pick transmitted DD grid if available (ground truth)
if isfield(W,'X')
    X_tx = W.X;
else
    X_tx = [];
end

% pick transmitted time-domain signal
if isfield(W,'tx_signal')
    tx = W.tx_signal(:);
elseif isfield(W,'tx')
    tx = W.tx(:);
else
    tx = [];
end

% pick received time-domain (prefer saved rx, else reconstruct approx from Xhat)
if isfield(R,'rx')
    rx = R.rx(:);
elseif isfield(R,'rx_signal')
    rx = R.rx_signal(:);
elseif isfield(R,'Y') && isfield(W,'params')
    % we might reconstruct but prefer rx; we'll try reconstructing from Xhat if available
    rx = [];
elseif isfield(R,'Xhat') && ~isempty(tx)
    % approximate rx by modulating Xhat through an approximate channel? skip — better to use Y
    rx = [];
else
    rx = [];
end

% if receiver_results contains Y (DD domain) use it directly
if isfield(R,'Y')
    Y_rx_dd = R.Y; % may be TF or DD depending on saved convention; we'll treat as DD from demodulator
else
    Y_rx_dd = [];
end

% If we don't have direct DD for tx or rx, use tx/rx time signals and demodulate to DD via SFFT
% define helper: time-domain to DD-domain demodulation (removes CP, reshape, SFFT)
function DD = td_to_dd(signal, Mlocal, Nlocal, Mcplocal)
    sig = signal(:);
    if Mcplocal > 0
        if numel(sig) > Mcplocal
            sig = sig(Mcplocal+1:end); % remove CP
        else
            warning('td_to_dd: signal shorter than Mcp; using full signal.');
        end
    end
    len = Mlocal * Nlocal;
    if numel(sig) < len
        sig = [sig; zeros(len - numel(sig),1)];
    elseif numel(sig) > len
        sig = sig(1:len);
    end
    % reshape - we must pick consistent ordering (N rows x M cols used in run_integration)
    Y_tf = reshape(sig, [Nlocal, Mlocal]);  % N x M
  % --- FIX: Match new otfs_demodulate (FFT Doppler, IFFT Delay) ---
DD = ifft(fft(Y_tf, [], 1), [], 2);
end

% Get Mcp if present
if isfield(params,'Mcp'), Mcp = params.Mcp; else Mcp = 0; end

% --- Build DD for transmitted signal ---
if ~isempty(X_tx)
    % We already have DD-grid X. But ensure orientation: our functions assume size N x M for Y (rows N, cols M)
    % X from Person A might be MxN or N x M; infer by matching dims
    [rX, cX] = size(X_tx);
    if rX == M && cX == N
        % X is M x N (delay bins x doppler?) convert to expected N x M by transpose
        Xtx_dd = X_tx.';
        warning('Assuming X in waveform_package is M x N; transposing to N x M for plotting.');
    elseif rX == N && cX == M
        Xtx_dd = X_tx;
    else
        % ambiguous: force reshape to N x M
        Xtx_dd = reshape(X_tx, [N, M]);
        warning('Reshaped X to [N M] convention for DD plotting.');
    end
else
    if ~isempty(tx)
        Xtx_dd = td_to_dd(tx, M, N, Mcp);
    else
        error('No transmitted DD grid (X) and no tx time-domain signal found in waveform_package.mat');
    end
end

% --- Build DD for received signal ---
if ~isempty(Y_rx_dd)
    % Received DD already available
    Ydd = Y_rx_dd;
else
    if ~isempty(rx)
        Ydd = td_to_dd(rx, M, N, Mcp);
    elseif isfield(R,'Xhat')
        % use estimated DD-domain Xhat from receiver results as proxy of received DD (after equalizer)
        Xhat = R.Xhat;
        % ensure shape
        [rXh, cXh] = size(Xhat);
        if rXh==M && cXh==N
            Ydd = Xhat.'; % transpose
        elseif rXh==N && cXh==M
            Ydd = Xhat;
        else
            Ydd = reshape(Xhat, [N, M]);
        end
        warning('Using Xhat as proxy for received DD domain representation (not raw received).');
    else
        error('No received signal available: please include rx or Y (DD) in receiver_results.mat or rx in run results.');
    end
end

% ---- Normalize for plotting (so magnitudes are comparable) ----
epsv = 1e-12;
mag_tx = abs(Xtx_dd);
mag_rx = abs(Ydd);

% Normalize to same peak (so we compare shapes rather than absolute gains)
mag_tx_n = mag_tx / (max(mag_tx(:))+epsv);
mag_rx_n = mag_rx / (max(mag_rx(:))+epsv);

% difference map
diff_map = mag_rx_n - mag_tx_n;

% ---- Basic diagnostics ----
fprintf('DD grids sizes used (rows x cols): tx: %d x %d, rx: %d x %d\n', size(Xtx_dd,1), size(Xtx_dd,2), size(Ydd,1), size(Ydd,2));

% Peak locations for TX and RX
[~, tx_idx] = max(mag_tx_n(:));
[tx_r, tx_c] = ind2sub(size(mag_tx_n), tx_idx);
[~, rx_idx] = max(mag_rx_n(:));
[rx_r, rx_c] = ind2sub(size(mag_rx_n), rx_idx);
fprintf('Strongest TX peak at (row=%d, col=%d). Strongest RX peak at (row=%d, col=%d)\n', tx_r, tx_c, rx_r, rx_c);

% ---- Plots ----
figure('Name','TX DD magnitude (normalized)','NumberTitle','off');
imagesc(mag_tx_n); axis xy; colorbar;
title('Transmitted: |X_{DD}| (normalized)');
xlabel('Delay bin (M dimension)');
ylabel('Doppler bin (N dimension)');

figure('Name','RX DD magnitude (normalized)','NumberTitle','off');
imagesc(mag_rx_n); axis xy; colorbar;
title('Received: |Y_{DD}| (normalized)');
xlabel('Delay bin (M dimension)');
ylabel('Doppler bin (N dimension)');

figure('Name','Difference (RX - TX) in DD (normalized)','NumberTitle','off');
imagesc(diff_map); axis xy; colorbar;
title('Normalized difference: |Y| - |X|');
xlabel('Delay bin (M dimension)');
ylabel('Doppler bin (N dimension)');

% Overlay peaks (if you want to mark multiple peaks)
% find top K peaks in RX
K = min(10, 50);
RXvals = mag_rx_n;
txvals = mag_tx_n;
% simple local maxima detection by sorting (non-ideal but sufficient)
[~, idxs] = sort(RXvals(:),'descend');
peaks = zeros(K,2);
for k=1:K
    [r,c] = ind2sub(size(RXvals), idxs(k));
    peaks(k,:) = [r,c];
end

figure('Name','RX DD with top-K peaks','NumberTitle','off');
imagesc(mag_rx_n); axis xy; colorbar; hold on;
plot(peaks(:,2), peaks(:,1), 'rx','MarkerSize',10,'LineWidth',1.5);
title(sprintf('Received DD magnitude with top-%d peaks',K));
xlabel('Delay bin (M dimension)'); ylabel('Doppler bin (N dimension)');
hold off;

% ---- Optional: show slices (delay slice at strongest Doppler and vice versa) ----
[~, dop_idx] = max(max(mag_rx_n,[],1)); % column index of strongest
[~, del_idx] = max(max(mag_rx_n,[],2)); % row index of strongest

figure('Name','Delay slice at strongest Doppler (TX vs RX)');
plot(1:size(mag_tx_n,1), mag_tx_n(:,dop_idx), 'b-','LineWidth',1.5); hold on;
plot(1:size(mag_rx_n,1), mag_rx_n(:,dop_idx), 'r--','LineWidth',1.5); hold off;
legend('TX','RX'); xlabel('Delay bin'); ylabel('Normalized magnitude');
title('Delay slice at strongest Doppler column');

figure('Name','Doppler slice at strongest Delay (TX vs RX)');
plot(1:size(mag_tx_n,2), mag_tx_n(del_idx,:), 'b-','LineWidth',1.5); hold on;
plot(1:size(mag_rx_n,2), mag_rx_n(del_idx,:), 'r--','LineWidth',1.5); hold off;
legend('TX','RX'); xlabel('Doppler bin'); ylabel('Normalized magnitude');
title('Doppler slice at strongest Delay row');

% ---- Save images to files ----
try
    saveas(figure(1),'tx_dd_magnitude.png');
    saveas(figure(2),'rx_dd_magnitude.png');
    saveas(figure(3),'dd_difference.png');
    saveas(figure(4),'rx_topk_peaks.png');
    saveas(figure(5),'delay_slice.png');
    saveas(figure(6),'doppler_slice.png');
    fprintf('Saved figures to current folder (png).\n');
catch
    fprintf('Could not save some figures (possibly due to figure numbering differences). Still plots shown.\n');
end

% ---- Print short guidance ----
fprintf('Interpretation tips:\n');
fprintf(' - If RX DD map is a shifted & scaled version of TX DD map: channel mainly applies delay/doppler shifts.\n');
fprintf(' - If RX shows spreading around TX impulses: fractional delays/dopplers or leakage / inter-Doppler interference.\n');
fprintf(' - If many spurious peaks appear: noise, multipath, or false detection (consider CFAR / thresholding).\n');

% --- end ---
