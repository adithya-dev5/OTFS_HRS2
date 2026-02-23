%% generate_ber_curve.m
% Generates a BER vs SNR graph for the corrected OTFS system.
% Uses "Perfect CSI" (Heff built from known targets) to verify performance.

clear; close all; clc;

% --- 1. SETUP & PATHS ---
if exist('src/receiver','dir'), addpath('src/receiver'); end
if exist('src/txr','dir'), addpath('src/txr'); end

% Simulation Parameters
M = 32; N = 32; 
Mcp = 0;           % Keep 0 for perfect alignment
Xp = 10;           % Strong Pilot
fs = 15.36e6;      % Sampling Rate
SNR_dB_Range = 0:2:20; % SNR Points to simulate
Max_Frames = 500;  % Frames per SNR point
Min_Errors = 300;  % Stop if enough errors are collected (for speed)

% --- 2. DEFINE STATIC CHANNEL (One Target) ---
% We use a static channel so we only build Heff ONCE (saves huge time)
params = struct('M',M, 'N',N, 'Mcp',Mcp, 'fs',fs, 'Lfft',2048);

% Target: 5 samples delay, 4 bins Doppler
d_samples = 5; 
k_bins = 4;
dt = 1/fs; 
df = fs/(M*N);

t1 = struct();
t1.gain = 0.8; 
t1.tau = d_samples * dt; 
t1.nu  = k_bins * df;
params.targets = t1;

fprintf('--- SYSTEM CONFIGURATION ---\n');
fprintf('Grid: %dx%d, Pilot: %d\n', M, N, Xp);
fprintf('Target: Delay=%d samples, Doppler=%d bins\n', d_samples, k_bins);

% --- 3. PRE-CALCULATE Heff (Perfect CSI) ---
fprintf('\nPre-building Effective Channel Matrix (Heff)... ');
% Create a dummy tx to initialize internal params if needed
dummy_tx = zeros(M*N, 1);
if exist('build_Heff','file')
    Heff = build_Heff(N, M, params, dummy_tx);
else
    error('build_Heff.m not found! Cannot equalize.');
end
fprintf('Done.\n');

% --- 4. MONTE CARLO LOOP ---
BER = zeros(size(SNR_dB_Range));

fprintf('\nStarting Simulation:\n');
fprintf('| SNR (dB) |  Frames  |  Errors  |   BER    |\n');
fprintf('|----------|----------|----------|----------|\n');

for i = 1:length(SNR_dB_Range)
    snr = SNR_dB_Range(i);
    
    total_errors = 0;
    total_bits = 0;
    frame_count = 0;
    
    while (total_errors < Min_Errors && frame_count < Max_Frames)
        frame_count = frame_count + 1;
        
        % A. GENERATE DATA
        % Pilot + Guard + Random QPSK
        [dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
        lp = 15; kp = 15;
        maskGuard = (abs(dd_L-lp)<=2) | (abs(dd_K-kp)<=2);
        maskData = ~maskGuard;
        
        DataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
        X = zeros(N, M);
        X(kp+1, lp+1) = Xp; % Pilot
        X(maskData) = DataSyms;
        
        % B. TRANSMIT (Modulate + Channel)
        tx = otfs_modulate(X, N, M);
        rx_clean = apply_channel(tx, params.targets, params);
        
        % C. ADD NOISE
        sig_pwr = mean(abs(rx_clean).^2);
        noise_pwr = sig_pwr / (10^(snr/10));
        noise = sqrt(noise_pwr/2) * (randn(size(rx_clean)) + 1j*randn(size(rx_clean)));
        rx_noisy = rx_clean + noise;
        
        % D. RECEIVE (Demodulate + Equalize)
        Y = otfs_demodulate(rx_noisy, Mcp, N, M);
        
        % MMSE Equalizer (using pre-built Heff)
        Xhat_vec = mmse_equalize(Y, Heff, noise_pwr);
        Xhat = reshape(Xhat_vec, N, M);
        
        % Normalize Xhat (Pilot Aided)
        pilot_est = Xhat(kp+1, lp+1);
        if abs(pilot_est) > 0.01
             Xhat = Xhat * (Xp / pilot_est);
        end
        
        % E. COUNT ERRORS
        % Extract Data Symbols
        rx_syms = Xhat(maskData);
        tx_syms = X(maskData);
        
        % Hard Decision (QPSK)
        rx_bits = qamdemod(rx_syms, 4, 'UnitAveragePower', true);
        tx_bits = qamdemod(tx_syms, 4, 'UnitAveragePower', true);
        
        % Bit Errors (2 bits per symbol)
        [~, ber_ratio] = biterr(tx_bits, rx_bits, 2); 
        errs = ber_ratio * (2 * numel(tx_bits));
        
        total_errors = total_errors + errs;
        total_bits = total_bits + (2 * numel(tx_bits));
    end
    
    BER(i) = total_errors / total_bits;
    fprintf('|   %4.1f   |   %4d   |   %5d  | %.2e |\n', snr, frame_count, total_errors, BER(i));
    
    % Early exit if BER is 0 (at high SNR)
    if total_errors == 0
        break; 
    end
end

% --- 5. PLOT ---
figure('Name', 'BER vs SNR', 'Color', 'w');
semilogy(SNR_dB_Range, BER, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('OTFS BER vs SNR');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-6 1]);
% --- 6. PLOTTING WITH UNITS ---
f = figure('Name', 'Annotated RVM', 'Color', 'w', 'Position', [100, 100, 900, 600]);

% Create Axes Vectors in PHYSICAL UNITS
% Shift the RVM matrix so the Pilot is at the center index.
pilot_idx_r = kp * Q + 1;
pilot_idx_c = lp * Q + 1;
RVM_shifted = circshift(RVM, [size(RVM,1)/2 - pilot_idx_r, size(RVM,2)/2 - pilot_idx_c]);

% Define Axes centered at 0
d_axis = ((1:size(RVM,2)) - size(RVM,2)/2) / Q * range_res; % Meters
v_axis = ((1:size(RVM,1)) - size(RVM,1)/2) / Q * vel_res;   % m/s

imagesc(d_axis, v_axis, RVM_shifted);
axis xy; colormap(jet); colorbar;
xlabel('Range (meters)');
ylabel('Velocity (m/s)');
title('Figure 5.2: Sensing Verification (Annotated)');

% --- FIX: MAKE GRID VISIBLE ---
grid on;
set(gca, 'Layer', 'top');       % Bring grid lines to the front
set(gca, 'GridColor', 'w');     % Make grid lines White
set(gca, 'GridAlpha', 0.4);     % Set transparency (0.4 is subtle but visible)
set(gca, 'LineWidth', 1.0);     % Make them slightly thicker
% ------------------------------
