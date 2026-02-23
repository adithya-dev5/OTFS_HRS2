% regenerate_waveform.m (Fixed for Visibility & Stability)
% Moves target away from Pilot so it appears clearly in RVM.

clear; clc; close all;
if exist('src/txr','dir'), addpath('src/txr'); end

% --- PARAMS ---
M = 32; N = 32; 
Mcp = 0;        % Keep 0 for clean alignment tests
fs = 1e6; 
Xp = 10;        % FIX: Increased from 1 to 10 so Pilot is stronger than Data
lp = 16; kp = 16; % Pilot at Center (0-based indices)
lTau = 2; kNu = 2;

% --- DATA GENERATION ---
rng(42);
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskGuard = (abs(dd_L-lp)<=lTau) | (abs(dd_K-kp)<=2*kNu);
maskData = ~maskGuard;

% Generate QPSK Data (Magnitude ~1)
% Note: Pilot (Mag 10) will be much stronger than these
DataSymbols = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);

% Build Grid
% Note: build_dd_grid expects 1-based indices, so we pass lp+1, kp+1
X = build_dd_grid(M, N, lp+1, kp+1, lTau, kNu, DataSymbols, Xp);

% Modulate
tx_signal = otfs_modulate(X, N, M);
if Mcp>0, tx_signal = [tx_signal(end-Mcp+1:end); tx_signal]; end

% --- TARGET DEFINITION ---
params = struct('M',M, 'N',N, 'Mcp',Mcp, 'fs',fs, 'Xp',Xp, 'lp',lp, 'kp',kp);
params.Lfft = max(1024, length(tx_signal)); % FIX: Define Lfft for channel accuracy

% Move Target to Bin (5, 8) relative to Pilot
% Resolution (High-Res Physics):
dt = 1/fs;      % 1 Sample = 1 Delay Bin
df = fs/(M*N);  % Standard Doppler Bin

t1 = struct();
t1.gain = 0.5;   % Weaker than pilot
t1.tau = 5 * dt; % Delay = 5 Samples (exactly 5 Bins in True OTFS)
t1.nu  = 8 * df; % Doppler = 8 Bins

% Add aliases for robustness
t1.delay = t1.tau; t1.doppler = t1.nu;

params.targets = [t1];

% --- SAVE ---
save('waveform_package.mat', 'X', 'tx_signal', 'params');
fprintf('SUCCESS: Waveform regenerated. Target at Relative Bin (5, 8).\n');