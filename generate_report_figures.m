% generate_report_figures_v2.m
% Generates Figures 4.1 to 4.7 for the End-Semester Report.
% SEPARATES "Sensing Physics" (Fractional) from "Comm Benchmark" (Integer).

clear; clc; close all;

% --- SETUP PATHS ---
if exist('src/receiver','dir'), addpath('src/receiver'); end
if exist('src/txr','dir'), addpath('src/txr'); end

% --- GLOBAL PARAMETERS ---
M = 32; N = 32; Mcp = 0; fs = 1e6; Xp = 1;
lp = 16; kp = 16; % Pilot Center (Indices: 17, 17)
lTau = 2; kNu = 2; % Guards

% --- DATA GENERATION ---
rng(42); 
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskPilot = (dd_L == lp) & (dd_K == kp);
maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );
maskData = ~(maskGuard | maskPilot);
% Use random QPSK data
DataSymbols = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);

% 1. Build Base Grid (X)
X = build_dd_grid(M, N, lp+1, kp+1, lTau, kNu, DataSymbols, Xp);
tx_sig = otfs_modulate(X, N, M);

% ---------------------------------------------------------
% FIG 4.1: Transmitted Single-Bin Pilot Structure
% ---------------------------------------------------------
fprintf('Generating Fig 4.1 (Tx Grid)...\n');
f1 = figure('Name','Fig 4.1','Visible','on');
imagesc(0:N-1, 0:M-1, abs(X)); axis xy; colorbar;
title('4.1 Transmitted Pilot Structure');
xlabel('Delay Bin'); ylabel('Doppler Bin');
caxis([0 1.2]);
saveas(f1, 'Fig4_1_TxGrid.png');

% =========================================================
% SCENARIO A: SENSING PHYSICS (Fractional Target)
% Used for Fig 4.2 (Raw Rx) and 4.3 (RVM)
% =========================================================
dt = 1/fs; df = fs/(M*N);
% Two Paths: LoS (Strong) + Fractional Target (Weak)
targets_frac = [];
targets_frac(1).gain = 1.0; targets_frac(1).tau = 0; targets_frac(1).nu = 0; 
targets_frac(2).gain = 0.5; targets_frac(2).tau = 5.5*dt; targets_frac(2).nu = 8.5*df; 
for i=1:2, targets_frac(i).delay=targets_frac(i).tau; targets_frac(i).doppler=targets_frac(i).nu; end
params_frac = struct('M',M,'N',N,'Mcp',Mcp,'fs',fs, 'targets', targets_frac);

% Channel Output
rx_frac = synth_channel_from_targets(tx_sig, params_frac);
Y_frac_raw = otfs_demodulate(rx_frac, Mcp, N, M);

% ROBUST ALIGNMENT: Force align to LoS (Pilot position)
% We expect peak at (kp+1, lp+1). Find peak in Y ONLY near that region to avoid locking to target.
pilot_region = Y_frac_raw(kp-2:kp+4, lp-2:lp+4); 
[~, idxLocal] = max(abs(pilot_region(:)));
[k_loc, l_loc] = ind2sub(size(pilot_region), idxLocal);
% Map back to global indices
k_peak = (kp-2) + k_loc - 1;
l_peak = (lp-2) + l_loc - 1;

% Calculate Shift (Targeting k=17, l=17)
shift_r = (kp+1) - k_peak;
shift_c = (lp+1) - l_peak;
Y_frac_aligned = circshift(Y_frac_raw, [shift_r, shift_c]);

% ---------------------------------------------------------
% FIG 4.2: Demodulated DD Map (Fractional/Aligned)
% ---------------------------------------------------------
fprintf('Generating Fig 4.2 (Aligned Rx - Fractional)...\n');
f2 = figure('Name','Fig 4.2','Visible','on');
imagesc(0:N-1, 0:M-1, abs(Y_frac_aligned)); axis xy; colorbar;
title('4.2 Aligned Rx DD Map (Fractional Channel)');
xlabel('Delay Bin'); ylabel('Doppler Bin');
saveas(f2, 'Fig4_2_RxGrid_Aligned.png');

% ---------------------------------------------------------
% FIG 4.3: High-Resolution RVM (Interpolated)
% ---------------------------------------------------------
fprintf('Generating Fig 4.3 (RVM)...\n');
% Dual-Domain Zero-Padding
Y_dual = fft2(Y_frac_aligned);
Y_dual_shifted = fftshift(Y_dual);
upsample = 4;
M_new = upsample * M; N_new = upsample * N;
Y_padded = zeros(M_new, N_new);
r_start = floor((M_new - M)/2) + 1; c_start = floor((N_new - N)/2) + 1;
Y_padded(r_start:r_start+M-1, c_start:c_start+N-1) = Y_dual_shifted;
Y_highres = ifft2(ifftshift(Y_padded));
RVM = abs(Y_highres); RVM = RVM / max(RVM(:));

f3 = figure('Name','Fig 4.3','Visible','on');
dop_axis = linspace(0, N-1, N_new); del_axis = linspace(0, M-1, M_new);
imagesc(dop_axis, del_axis, RVM); axis xy; colormap(jet); colorbar;
title('4.3 High-Resolution RVM');
xlabel('Doppler Bin'); ylabel('Delay Bin');
xlim([0 N-1]); ylim([0 M-1]);
saveas(f3, 'Fig4_3_RVM.png');

% =========================================================
% SCENARIO B: COMM BENCHMARK (Integer Channel)
% Used for Fig 4.4 - 4.7 (Clean Constellation, Low BER)
% =========================================================
% Integer Target (Aligned) to ensure clean equalization
targets_int = targets_frac;
targets_int(2).tau = 5*dt; targets_int(2).nu = 8*df; % Integer bins
for i=1:2, targets_int(i).delay=targets_int(i).tau; targets_int(i).doppler=targets_int(i).nu; end
params_int = params_frac; params_int.targets = targets_int;

rx_int = synth_channel_from_targets(tx_sig, params_int);
Y_int = otfs_demodulate(rx_int, Mcp, N, M);
Y_int_aligned = circshift(Y_int, [shift_r, shift_c]); % Use same shift

% ---------------------------------------------------------
% FIG 4.4: Demodulated DD Map (Integer)
% ---------------------------------------------------------
fprintf('Generating Fig 4.4 (Integer Channel)...\n');
f4 = figure('Name','Fig 4.4','Visible','on');
imagesc(0:N-1, 0:M-1, abs(Y_int_aligned)); axis xy; colorbar;
title('4.4 Integer-Channel DD Map');
xlabel('Delay Bin'); ylabel('Doppler Bin');
saveas(f4, 'Fig4_4_RxGrid_Integer.png');

% ---------------------------------------------------------
% FIG 4.5: Equalized DD Grid (MMSE on Integer Channel)
% ---------------------------------------------------------
fprintf('Generating Fig 4.5 (Equalized Grid)...\n');
Heff = build_Heff(N, M, params_int, tx_sig); % Perfect CSI
Xhat_vec = mmse_equalize(Y_int, Heff, 1e-4); % Use RAW Y_int
Xhat = reshape(Xhat_vec, N, M);

% Normalize (Pilot Peak)
pilot_est = Xhat(kp+1, lp+1);
if abs(pilot_est) > 1e-9, Xhat = Xhat * (Xp / pilot_est); end

f5 = figure('Name','Fig 4.5','Visible','on');
imagesc(0:N-1, 0:M-1, abs(Xhat)); axis xy; colorbar;
title('4.5 Equalized DD Grid');
xlabel('Delay Bin'); ylabel('Doppler Bin');
caxis([0 1.2]);
saveas(f5, 'Fig4_5_EqualizedGrid.png');

% ---------------------------------------------------------
% FIG 4.6: Detection Probability vs SNR (Monte Carlo)
% ---------------------------------------------------------
fprintf('Generating Fig 4.6 (Pd vs SNR)...\n');
snr_range = -10:2:10; pd_log = zeros(size(snr_range));
num_trials = 50;
% Expected Peak in Aligned Grid (Integer Target: 5, 8 relative to Pilot 16,16)
exp_k = kp + 1 + 8; exp_l = lp + 1 + 5; 

for i = 1:length(snr_range)
    n_pwr = mean(abs(rx_int).^2) / (10^(snr_range(i)/10));
    detects = 0;
    for t = 1:num_trials
        noise = sqrt(n_pwr/2)*(randn(size(rx_int))+1j*randn(size(rx_int)));
        Y_n = otfs_demodulate(rx_int + noise, Mcp, N, M);
        % Align (Fast)
        [~, iY] = max(abs(Y_n(:))); [ky, ly] = ind2sub(size(Y_n), iY);
        Y_a = circshift(Y_n, [(kp+1)-ky, (lp+1)-ly]);
        % Mask
        Y_search = abs(Y_a); Y_search(maskGuard|maskPilot) = 0;
        [~, iT] = max(Y_search(:)); [kT, lT] = ind2sub(size(Y_search), iT);
        if abs(kT - exp_k) <= 1 && abs(lT - exp_l) <= 1, detects = detects + 1; end
    end
    pd_log(i) = detects / num_trials;
end

f6 = figure('Name','Fig 4.6','Visible','on');
plot(snr_range, pd_log, 'r-o', 'LineWidth', 2); grid on;
title('4.6 Detection Probability vs SNR');
xlabel('SNR (dB)'); ylabel('Pd'); ylim([-0.1 1.1]);
saveas(f6, 'Fig4_6_Pd_vs_SNR.png');

% ---------------------------------------------------------
% FIG 4.7: Recovered Data Symbols (Constellation)
% ---------------------------------------------------------
fprintf('Generating Fig 4.7 (Constellation)...\n');
% Use Xhat from Fig 4.5 (Integer Channel = Clean Constellation)
syms_rx = Xhat(maskData);
syms_tx = X(maskData);

f7 = figure('Name','Fig 4.7','Visible','on');
plot(real(syms_rx), imag(syms_rx), 'b.', 'MarkerSize', 8); hold on;
plot(real(syms_tx), imag(syms_tx), 'rx', 'LineWidth', 2, 'MarkerSize', 10);
grid on; axis equal;
title('4.7 Recovered Data Symbols');
legend('Received (Equalized)','Transmitted');
saveas(f7, 'Fig4_7_Constellation.png');

fprintf('DONE. All 7 Figures Saved.\n');
% Optional: Visualize the Effective Channel Matrix (Heff)
figure('Name', 'Extra: Heff Matrix Structure');
imagesc(abs(Heff)); 
axis xy; axis square; colorbar;
title('Effective Channel Matrix H_{eff} (Magnitude)');
xlabel('Transmit Index (n,m)');
ylabel('Receive Index (l,k)');
fprintf('Heff visualization generated.\n');

% ... [Append this to the end of your existing script] ...

% ---------------------------------------------------------
% FIG 4.8: Sensing Performance (Detection Probability vs SNR)
% ---------------------------------------------------------
fprintf('Generating Fig 4.8 (Sensing Pd vs SNR)...\n');

% Parameters for Sensing Benchmark
snr_range = -10:2:15; 
pd_log = zeros(size(snr_range));
num_trials = 100; % 100 Monte Carlo runs per SNR point

% Target Location (Integer Aligned for fairness)
% Pilot at (16, 16). Target at (21, 24) -> (+5 Delay, +8 Doppler)
exp_k = kp + 1 + 8; 
exp_l = lp + 1 + 5; 

% Calculate signal power from the integer reference signal
sig_pwr = mean(abs(rx_int).^2);

for i = 1:length(snr_range)
    % Calculate Noise Power for current SNR
    n_pwr = sig_pwr / (10^(snr_range(i)/10));
    detects = 0;
    
    for t = 1:num_trials
        % Add Noise
        noise = sqrt(n_pwr/2) * (randn(size(rx_int)) + 1j*randn(size(rx_int)));
        Y_n = otfs_demodulate(rx_int + noise, Mcp, N, M);
        
        % Alignment (Using your V2/V6 robust logic)
        pilot_reg = Y_n(kp-2:kp+4, lp-2:lp+4);
        [~, idxL] = max(abs(pilot_reg(:))); 
        [kl, ll] = ind2sub(size(pilot_reg), idxL);
        k_pk = (kp-2)+kl-1; 
        l_pk = (lp-2)+ll-1;
        Y_a = circshift(Y_n, [(kp+1)-k_pk, (lp+1)-l_pk]);
        
        % Mask Pilot & Guard (Don't detect the pilot as a target!)
        Y_search = abs(Y_a); 
        Y_search(maskGuard | maskPilot) = 0;
        
        % Peak Detection
        [~, iT] = max(Y_search(:)); 
        [kT, lT] = ind2sub(size(Y_search), iT);
        
        % Hit? (Allow +/- 1 bin tolerance for noise jitter)
        if abs(kT - exp_k) <= 1 && abs(lT - exp_l) <= 1
            detects = detects + 1;
        end
    end
    pd_log(i) = detects / num_trials;
    % fprintf('SNR %d dB: Pd = %.2f\n', snr_range(i), pd_log(i)); 
end

% Plot
f8 = figure('Name','Fig 4.8','Visible','on');
plot(snr_range, pd_log, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r'); 
grid on;
title('4.8 Detection Probability (Pd) vs SNR');
xlabel('SNR (dB)'); ylabel('Detection Probability (Pd)');
ylim([-0.05 1.05]);
yline(0.9, 'k--', '90% Reliability'); % Reference line
legend('OTFS Sensing', '90% Threshold', 'Location', 'SouthEast');
saveas(f8, 'Fig4_8_Pd_vs_SNR.png');

fprintf('DONE. Fig 4.8 Generated.\n');