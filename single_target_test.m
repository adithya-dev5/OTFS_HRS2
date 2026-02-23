% single_target_test.m
% Sanity test: 1 DD-bin -> modulate -> apply_channel -> demodulate -> inspect Y.

clear; close all; clc;

if exist('src/txr','dir')
    addpath('src/txr');
end

% Load waveform_package to get params and convenience
wpfile = 'waveform_package.mat';
if ~isfile(wpfile), error('Put waveform_package.mat in project root before running this test.'); end
S = load(wpfile);
params = S.params;

% Ensure required params exist or set reasonable defaults:
if ~isfield(params,'N') || ~isfield(params,'M')
    error('params.N and params.M must exist in waveform_package.mat');
end
N = params.N; M = params.M;

% Ensure Delta_f and fs present (required to convert indices -> tau/nu)
if ~isfield(params,'Delta_f')
    if isfield(params,'T')
        params.Delta_f = 1/params.T;
    else
        % choose a default subcarrier spacing if not provided (adjust if needed)
        params.Delta_f = 15000; % 15 kHz example (set to your system value)
        warning('params.Delta_f not found. Using default Delta_f = %g Hz', params.Delta_f);
    end
end
if ~isfield(params,'fs')
    % sampling freq — choose M * Delta_f as typical OTFS sampling grid (one TF frame)
    params.fs = params.M * params.Delta_f;
    warning('params.fs not found. Using params.fs = M*Delta_f = %g Hz', params.fs);
end

% Choose a target index in DD (li = delay index along columns, ki = doppler index along rows)
% NOTE: your X convention: otfs_modulate assumed X is size [N, M] (rows = N (doppler), cols = M (delay)).
ki_test = 5;   % Doppler index (1..N)
li_test = 3;   % Delay index (1..M)
gain = 1+0j;   % complex gain

% Build a single-DD bin transmit grid X (N x M)
X_single = zeros(N, M);
% place a bin at (row=ki_test, col=li_test) — check your indexing/orientation convention
X_single(ki_test, li_test) = 1;

% Produce time-domain transmit signal (uses your otfs_modulate)
tx = otfs_modulate(X_single, N, M);

% Prepare params.targets with tau, nu, gain required by apply_channel
% Conversion formulas:
%   tau = li / (M * Delta_f)
%   nu  = ki * Delta_f / N
tau = (li_test-1) / (M * params.Delta_f);   % subtract 1 if your indices are 1-based; adjust if your 'li' was 0-based
nu  = (ki_test-1) * params.Delta_f / N;     % same adjustment for 1-based index
% If your li,ki are already 0-based, remove the -1.

% Create targets struct with required fields
tg = struct();
tg(1).tau  = tau;     % seconds
tg(1).nu   = nu;      % Hz
tg(1).gain = gain;

% ensure Lfft and fs defined for apply_channel
if ~isfield(params,'Lfft') || isempty(params.Lfft)
    params.Lfft = max(2^nextpow2(4*length(tx)), length(tx));
end

% Call apply_channel (preferred physical-channel)
rx = apply_channel(tx, tg, params);   % this expects tg with fields tau, nu, gain

% Add CP if your pipeline expects it (here we assume run_integration handled CP logic)
if ~isfield(params,'Mcp'), params.Mcp = 0; end
if params.Mcp>0
    rx_cp = [rx(end-params.Mcp+1:end); rx(:)];
else
    rx_cp = rx(:);
end

% Demodulate using your corrected otfs_demodulate (this returns Y in Delay x Doppler if yours transposes)
Y_test = otfs_demodulate(rx_cp, params.Mcp, N, M);

% Visualize
figure; imagesc(abs(X_single)); axis xy; colorbar; title('TX single DD bin |X| (grid)');
figure; imagesc(abs(Y_test)); axis xy; colorbar; title('RX demodulated |Y| — single-target test');

% Diagnostics: find peak in Y_test
[~, idxY] = max(abs(Y_test(:)));
[ry, cy] = ind2sub(size(Y_test), idxY);
fprintf('Peak in Y at (row=%d, col=%d)  -- expected (ki=%d, li=%d)\n', ry, cy, ki_test, li_test);

% Compute pilot/data energy match (if desired)
frac_in_tx_bin = sum(abs(Y_test(ki_test, li_test))^2) / (sum(abs(Y_test(:)).^2) + eps);
fprintf('Energy fraction at expected bin = %.4f\n', frac_in_tx_bin);

% Save Y_test for inspection
save('single_target_Y_test.mat','Y_test','X_single','tg','params');

