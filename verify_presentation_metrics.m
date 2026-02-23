%% verify_presentation_metrics.m
% Generates the exact numbers needed for Slides 17, 18, and 19.
% FIXED: Guard mask aligned with build_dd_grid logic (2*kNu).

clear; close all; clc;
if exist('src/txr','dir'), addpath('src/txr'); end
if exist('src/receiver','dir'), addpath('src/receiver'); end

% --- CONFIG ---
M = 32; N = 32; Mcp = 32;
fs = 15.36e6; Xp = 10;
lp = 15; kp = 15;
lTau = 2; kNu = 2; % Define Guard Params explicitly

params = struct('M',M, 'N',N, 'Mcp',Mcp, 'fs',fs, 'Lfft',2048);

% --- PART 1: VALIDATION METRICS (Slide 19) ---
fprintf('--- Generating Data for Slide 19 (Validation) ---\n');

% 1. Build Heff
params.targets = struct('tau',0, 'nu',0, 'gain',1); % Static Channel
Heff = build_Heff(N, M, params, []);

% 2. Noiseless Transmission
rng(42);
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);

% --- FIX: Match build_dd_grid logic (Doppler guard is 2*kNu) ---
maskGuard = (abs(dd_L-lp)<=lTau) | (abs(dd_K-kp)<=2*kNu); 
maskPilot = (dd_L==lp) & (dd_K==kp);
maskData = ~(maskGuard | maskPilot);

% Generate exactly the right number of symbols
DataSyms = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);

X = build_dd_grid(M, N, lp, kp, lTau, kNu, DataSyms, Xp);
tx = otfs_modulate(X, N, M);
if Mcp>0, tx = [tx(end-Mcp+1:end); tx]; end

% 3. Noiseless Reception
rx = apply_channel(tx, params.targets, params);
Y = otfs_demodulate(rx, Mcp, N, M);

% 4. Perfect Equalization
x_vec = mmse_equalize(Y, Heff, 0); 
Xhat = reshape(x_vec, N, M);

% 5. Calculate MSE & SER
err = X - Xhat;
mse_val = mean(abs(err(:)).^2);
ser_val = sum(abs(err(maskData)) > 1e-3);

fprintf('Noiseless MSE: %.2e\n', mse_val);
fprintf('Noiseless SER: %d\n', ser_val);
fprintf('-----------------------------------------------\n');


% --- PART 2: SNR SWEEP (Slide 17) ---
fprintf('\n--- Generating Data for Slide 17 (SNR Sweep) ---\n');
SNR_list = [0, 5, 10, 15, 20];
% High-Fidelity Fractional Target
params.targets = struct('gain', 0.8, 'tau', 5.3/fs, 'nu', 4.2*(fs/(M*N))); 

fprintf('SNR (dB) |   BER    | Peak Detection Status\n');
fprintf('---------|----------|----------------------\n');

sig_pwr = mean(abs(rx).^2);

for snr = SNR_list
    n_pwr = sig_pwr / (10^(snr/10));
    
    % Monte Carlo (Short run for speed)
    errs = 0; bits = 0;
    for f = 1:50 
        noise = sqrt(n_pwr/2) * (randn(size(rx)) + 1j*randn(size(rx)));
        y_noisy = otfs_demodulate(rx + noise, Mcp, N, M);
        x_eq = reshape(mmse_equalize(y_noisy, Heff, n_pwr), N, M);
        
        % Normalize
        pe = x_eq(kp+1, lp+1); 
        if abs(pe)>0.1, x_eq = x_eq * (Xp/pe); end
        
        % BER
        tx_s = X(maskData); rx_s = x_eq(maskData);
        [~, r] = biterr(qamdemod(tx_s,4), qamdemod(rx_s,4), 2);
        errs = errs + r*2*numel(tx_s);
        bits = bits + 2*numel(tx_s);
        if errs > 200, break; end
    end
    
    ber = errs/bits;
    status = 'Unclear';
    if ber < 1e-4, status = 'Perfect';
    elseif ber < 1e-2, status = 'Reliable'; 
    end
    
    fprintf('   %2d    | %.2e | %s\n', snr, ber, status);
end