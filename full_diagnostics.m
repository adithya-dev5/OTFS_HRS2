% generate_perfect_plots.m
% Generates "Textbook" OTFS results by snapping targets to Integer Grid.
% Removes Fractional Leakage to show the theoretical performance limit.

clear; clc; close all;

% --- SETUP ---
if exist('src/receiver','dir'), addpath('src/receiver'); end
if exist('src/txr','dir'), addpath('src/txr'); end

if ~isfile('waveform_package.mat'), error('waveform_package.mat missing.'); end
load('waveform_package.mat'); 

% --- STEP 1: SNAP TARGETS TO INTEGER GRID ---
% Calculate Resolution
doppler_res = params.fs / (params.M * params.N); % Hz per bin
delay_res   = 1 / params.fs;                     % Sec per bin

fprintf('Original Target Doppler: %.2f Hz\n', params.targets(1).nu);

% Modify Param to be EXACTLY Integer Bin 1 (removes leakage)
% We shift the physics to match the math.
params.targets(1).nu = doppler_res * 1.0;   % Exactly Bin 1
params.targets(1).doppler = params.targets(1).nu; % Alias
params.targets(1).tau = delay_res * 2.0;    % Exactly Bin 2
params.targets(1).delay = params.targets(1).tau;  % Alias

fprintf('Adjusted Target Doppler: %.4f Hz (Exact Integer Bin 1)\n', params.targets(1).nu);
fprintf('This removes fractional leakage for the "Perfect" plot.\n');

% Force Consistent Channel Model
channel_func = @(tx_in, p) synth_channel_from_targets(tx_in, p);


% --- STEP 2: REBUILD Heff (With New Target) ---
fprintf('\nBuilding Heff (Aligned)... ');
MN = params.M * params.N;
Heff = zeros(MN, MN);

for k = 1:MN
    x_imp_dd = zeros(MN,1); x_imp_dd(k) = 1;
    X_imp = reshape(x_imp_dd, params.N, params.M);
    tx_imp = otfs_modulate(X_imp, params.N, params.M);
    if params.Mcp > 0, tx_imp = [tx_imp(end-params.Mcp+1:end); tx_imp]; end
    
    % Use NEW params with aligned target
    rx_imp = channel_func(tx_imp, params);
    
    Y_imp = otfs_demodulate(rx_imp, params.Mcp, params.N, params.M);
    Heff(:, k) = Y_imp(:);
    if mod(k, 250)==0, fprintf('%.0f%%.. ', k/MN*100); end
end
fprintf('Done.\n');


% --- STEP 3: BER SIMULATION ---
fprintf('\nSimulating BER vs SNR (Integer Aligned)...\n');
snr_range = 0:4:24; % 0 to 24 dB
ber_log = zeros(size(snr_range));

% Prepare Transmit Signal (With New Target Channel)
tx_fresh = otfs_modulate(X, params.N, params.M);
if params.Mcp > 0, tx_fresh = [tx_fresh(end-params.Mcp+1:end); tx_fresh]; end
rx_clean = channel_func(tx_fresh, params); % Noise-free Rx
sig_pwr = mean(abs(rx_clean).^2);

% Masks
[~, idxPeak] = max(abs(X(:))); 
[kp, lp] = ind2sub(size(X), idxPeak);
maskAll = abs(X) > 0;
maskData = maskAll; maskData(kp, lp) = false;
num_bits_frame = nnz(maskData) * 2;

% Settings
MIN_ERRORS = 30; MAX_FRAMES = 200;

for i = 1:length(snr_range)
    snr_val = snr_range(i);
    fprintf('  SNR: %2d dB ... ', snr_val);
    
    n_pwr = sig_pwr / (10^(snr_val/10));
    total_errs = 0; total_bits = 0; frames = 0;
    
    while (total_errs < MIN_ERRORS && frames < MAX_FRAMES)
        frames = frames + 1;
        
        % Noise
        noise = sqrt(n_pwr/2) * (randn(size(rx_clean)) + 1j*randn(size(rx_clean)));
        rx_loop = rx_clean + noise;
        
        % Rx
        Y_raw = otfs_demodulate(rx_loop, params.Mcp, params.N, params.M);
        
        % Equalize (No manual shift needed, Heff handles integer shift perfectly)
        Xhat_vec = mmse_equalize(Y_raw, Heff, n_pwr);
        Xhat = reshape(Xhat_vec, params.N, params.M);
        
        % Normalize
        pilot_est = Xhat(kp, lp);
        if abs(pilot_est) > 1e-9, Xhat = Xhat * (X(kp,lp) / pilot_est); end
        
        % Count
        rx_bits = [real(Xhat(maskData))>0, imag(Xhat(maskData))>0]; 
        tx_bits = [real(X(maskData))>0, imag(X(maskData))>0];
        errs = sum(rx_bits(:) ~= tx_bits(:));
        
        total_errs = total_errs + errs;
        total_bits = total_bits + numel(tx_bits);
        
        % Early exit for high SNR (if 0 errors after 20 frames, assume 0)
        if snr_val >= 16 && frames >= 20 && total_errs == 0, break; end
    end
    
    if total_bits > 0, ber_log(i) = total_errs / total_bits; else, ber_log(i) = 0; end
    fprintf('BER = %.4e (Errs: %d)\n', ber_log(i), total_errs);
end

% --- PLOT ---
figure('Position',[100 100 600 500]);
semilogy(snr_range, ber_log, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('OTFS Performance (Integer-Aligned)');
ylim([1e-6 1]);
saveas(gcf, 'Perfect_BER_Curve.png');
fprintf('\nSuccess! Saved "Perfect_BER_Curve.png".\n');