%% otfs_isac_upgraded_figures.m
% Upgraded OTFS Transmitter + Receiver with separate figures for:
%  - Estimated DD grid (Y) magnitude
%  - Transmit signal magnitude
% Features:
%  - frequency-domain fractional delay (phase ramp)
%  - exact Doppler multiplication in time domain
%  - vectorized high-resolution RVM via zero-padded 2D FFT
%  - coarse per-bin equalization (placeholder for Heff/MMSE)
% Author: ChatGPT (adapted to user's project)

clear; close all; clc;

%% ---------------- Simulation Parameters ----------------
M = 32;           % # delay bins (columns)
N = 32;           % # Doppler bins (rows)
B = 20e6;         % bandwidth (Hz)
fc = 77e9;        % carrier frequency (Hz)
Delta_f = B/M;    % subcarrier spacing
T = 1/Delta_f;    % symbol duration
frameT = N*T;     % frame duration
Ts = 1/(M*Delta_f); % sample interval (1/B)
fs = 1/Ts;        % sampling frequency
MN = M*N;

% Pilot location (0-based theory -> +1 for MATLAB)
lp = floor(M/2); kp = floor(N/2);
lp_m = lp + 1; kp_m = kp + 1;

% Guard half-widths
lTau = 2; kNu = 2;

% Pilot amplitude
Xp = 10;

% Targets: [Range_m, Velocity_mps, complex_gain]
targets = [15.5, 40, 1+0.2j;
           20.0, 10.4, 0.8+0.1j;
           42.0,-55, 0.9+0.05j];
Ltargets = size(targets,1);

c = 3e8;

% Convert physical params -> continuous delay/doppler bins
delay_bin_seconds = 1/(M*Delta_f); % = Ts
doppler_bin_hz = 1/(N*T);

tau_cont = 2*targets(:,1)/c;                % round-trip delay (s)
nu_hz   = 2*targets(:,2)*fc/c;              % Doppler freq (Hz)

li_cont = tau_cont / delay_bin_seconds;     % delay index (float)
ki_cont = nu_hz / doppler_bin_hz;           % doppler index (float)

li_int = floor(li_cont); li_frac = li_cont - li_int;
ki_int = floor(ki_cont); ki_frac = ki_cont - ki_int;

% SNR for pilot (linear)
SNRp_dB = 20;
SNRp = 10^(SNRp_dB/10);

%% ---------------- Build Delay-Doppler Grid (N x M) ----------------
X = zeros(N, M); % rows: Doppler k, cols: Delay l

% Place pilot
X(kp_m, lp_m) = Xp;

% Guard (explicit)
for k = (kp_m - 2*kNu):(kp_m + 2*kNu)
    for l = (lp_m - lTau):(lp_m + lTau)
        if k>=1 && k<=N && l>=1 && l<=M
            if ~(k==kp_m && l==lp_m)
                X(k,l) = 0;
            end
        end
    end
end

% Fill remaining grid with QPSK data
data_mask = (X==0);
num_data = sum(data_mask(:));
rng(0);
data_bits = randi([0 1], num_data*2,1);
data_sym = (1-2*data_bits(1:2:end)) + 1j*(1-2*data_bits(2:2:end));
X(data_mask) = data_sym;

%% ---------------- OTFS Modulator (approx ISFFT via column-wise IFFT) ----------------
S = ifft(X, [], 1);        % N x M
tx_signal = reshape(S, [], 1); % column-major -> length N*M

% Add cyclic prefix
Mcp = 16; % CP length
tx_cp = [tx_signal(end-Mcp+1:end); tx_signal];
len_cp = length(tx_cp);

% normalize tx energy (helpful for consistent SNR)
tx_cp = tx_cp / std(tx_cp(:));
tx = tx_signal / std(tx_signal(:)); % normalized tx (no CP) for later plotting/use

%% ---------------- Channel: frequency-domain fractional delay + Doppler ----------------
tx = tx_signal; % base signal (no CP for channel ops)
Ltx = length(tx);

% FFT length for fractional delay (zero-pad)
Lfft = 2^nextpow2(2*Ltx);
fvec = (0:(Lfft-1)).' * (fs / Lfft); % Hz
tvec = (0:Ltx-1).' * Ts; % seconds

rx_nocp_sum = zeros(Ltx,1);

for p=1:Ltargets
    hi = targets(p,3);
    tau = tau_cont(p);       % seconds
    nu  = nu_hz(p);          % Hz
    
    TXF = fft(tx, Lfft);
    Hdelay = exp(-1j*2*pi * fvec * tau);    % freq-domain delay phase ramp
    tx_delay = ifft(TXF .* Hdelay, Lfft);
    tx_delay = tx_delay(1:Ltx);             % crop to original length
    
    tx_delay_doppler = tx_delay .* exp(1j*2*pi*nu*tvec); % apply Doppler in time domain
    
    path_sig = hi * tx_delay_doppler;
    rx_nocp_sum = rx_nocp_sum + path_sig;
end

% attach CP
rx_with_cp = [rx_nocp_sum(end-Mcp+1:end); rx_nocp_sum];

% AWGN based on pilot energy
pilot_energy = mean(abs(tx_cp).^2);
noise_var = pilot_energy / SNRp;
rng(1);
rx = rx_with_cp + sqrt(noise_var/2)*(randn(size(rx_with_cp)) + 1j*randn(size(rx_with_cp)));

%% ---------------- Receiver: remove CP, reshape, FFT (SFFT approx) ----------------
rx_nocp = rx(Mcp+1:end); % remove CP
Y_time = reshape(rx_nocp, N, M);  % N x M
Y = fft(Y_time, [], 1);           % N x M -> delay-Doppler domain estimate

% --- Figure: Magnitude of Estimated Delay–Doppler Grid (separate figure) ---
figure('Name','Estimated DD Grid (Y)');
imagesc(abs(Y)); axis xy;
xlabel('Delay bin (l)'); ylabel('Doppler bin (k)');
title('Magnitude of Estimated DD Grid (Y)');
colorbar;

% --- Figure: Transmit-Signal Magnitude (separate figure) ---
figure('Name','Transmit Signal Magnitude');
plot(abs(tx),'LineWidth',1);
xlabel('Sample index');
ylabel('Magnitude');
title('Transmit Signal Magnitude');
grid on;

%% ---------------- High-Resolution Range-Velocity Map (via zero-padded 2D FFT) ----------------
eps_theta = 100;  % finer in delay
eps_phi   = 100;  % finer in doppler

Nfine = N * eps_phi;
Mfine = M * eps_theta;

% 2D FFT interpolation via zero-padding
RVM_complex = fft2(Y, Nfine, Mfine);
RVM = abs(fftshift(RVM_complex));

% find top Ltargets peaks
RVM_copy = RVM;
est_fine_indices = zeros(Ltargets,2); % [col (range index), row (vel index)]
for t = 1:Ltargets
    [~, idx] = max(RVM_copy(:));
    [r_idx, c_idx] = ind2sub(size(RVM_copy), idx);
    est_fine_indices(t,:) = [c_idx, r_idx];
    RVM_copy(r_idx, c_idx) = 0;
end

% center indices for mapping
center_row = Nfine/2 + 1;
center_col = Mfine/2 + 1;

li_hat = zeros(Ltargets,1);
ki_hat = zeros(Ltargets,1);
for t=1:Ltargets
    c_idx = est_fine_indices(t,1);
    r_idx = est_fine_indices(t,2);
    rel_col = c_idx - center_col;
    rel_row = r_idx - center_row;
    li_coarse = mod(rel_col/eps_theta, M);
    ki_coarse = mod(rel_row/eps_phi, N);
    li_hat(t) = li_coarse;
    ki_hat(t) = ki_coarse;
end

% Convert to physical Range and Velocity using report formulas:
Ri_hat = (li_hat) * c ./ (2 * M * Delta_f);   % meters
Vi_hat = (ki_hat) * c ./ (2 * fc * N * T);    % m/s

fprintf('\nEstimated Ranges (m) from high-res RVM (unordered):\n');
disp(Ri_hat.');
fprintf('Estimated Velocities (m/s) from high-res RVM (unordered):\n');
disp(Vi_hat.');

% Plot RVM with markers
figure('Name','High-res Range-Velocity Map');
imagesc(linspace(-M/2, M/2, Mfine), linspace(-N/2, N/2, Nfine), RVM); axis xy;
xlabel('Fine-range index (signed)'); ylabel('Fine-velocity index (signed)');
title('High-resolution RVM (zero-padded 2D FFT)');
colorbar;
hold on;
for t=1:Ltargets
    plot(est_fine_indices(t,1), est_fine_indices(t,2), 'wx', 'MarkerSize', 10, 'LineWidth', 2);
end
hold off;

%% ---------------- Simple MMSE-like Equalization (coarse placeholder) ----------------
[~, idxYmax] = max(abs(Y(:)));
[kmax, lmax] = ind2sub(size(Y), idxYmax);
H_est_pilot = Y(kmax, lmax) / Xp;

H_est = H_est_pilot * ones(size(Y)); % naive constant channel
Xhat = Y ./ H_est;

% Recover data symbols (naive)
recov_symbols = Xhat(data_mask);
orig_symbols = X(data_mask);

% Crude symbol error rate (naive)
ser = mean(sign(real(recov_symbols)) ~= sign(real(orig_symbols)) | sign(imag(recov_symbols)) ~= sign(imag(orig_symbols)));
fprintf('\nCrude symbol-error-rate (naive equalizer): %.4f\n', ser);

%% ---------------- Final Transmit signal plot (if you want a bigger separate window) ----------------
% Already plotted above. End of script.


%% ---------------- Extension: Heff construction + MMSE, BER & RMSE vs SNR ----------------
% This block constructs Heff empirically, performs MMSE equalization,
% and evaluates BER and RMSE vs SNR.

% Parameters (tweakable)
SNRdB_list = 0:5:30;      % SNR range to sweep (dB)
nSNR = length(SNRdB_list);
nTrials = 20;             % Monte Carlo trials per SNR (reduce to speed up)

% Reuse variables from above script:
% M, N, Delta_f, T, Ts, fs, Mcp, tx_signal, tx_cp, tx, X, data_mask,
% targets, tau_cont, nu_hz, Lfft, fvec, tvec

MN = M * N;
fprintf('\nBuilding empirical Heff matrix (this may take a while, MN = %d)...\n', MN);

% Preallocate Heff: size (MN_rx x MN_tx) where MN_rx = MN (vectorized Y), MN_tx = MN (vectorized X)
MN_rx = MN; MN_tx = MN;
Heff = zeros(MN_rx, MN_tx);

% We'll need functions used earlier: mapping Xgrid -> tx_signal and channel & receiver ops.
% For speed, precompute fft basis lengths etc used earlier.

% Helper: function that maps an input DD-vector x_vec (length MN_tx) to y_vec (length MN_rx)
function y_vec = txch_rx_map(x_vec, N, M, Mcp, Lfft, fvec, tvec, targets, tx_norm_flag)
    % x_vec is vectorized X (size N*M) in column-major order (same convention used before)
    Xgrid = reshape(x_vec, N, M);
    % OTFS modulator (ISFFT approx): ifft along rows
    S_local = ifft(Xgrid, [], 1);
    tx_local = reshape(S_local, [], 1);
    % Normalize to same scale as earlier tx (optional)
    if tx_norm_flag
        tx_local = tx_local / std(tx_local(:));
    end
    Ltx_local = length(tx_local);
    % Frequency-domain delay operator uses Lfft that's >= 2*Ltx_local (we pass same Lfft)
    % compute TXF for local (zero-pad to Lfft)
    TXF_local = fft(tx_local, Lfft);
    % build rx sum
    rx_nocp_sum_local = zeros(Ltx_local,1);
    for pp = 1:size(targets,1)
        hi = targets(pp,3);
        tau = 2*targets(pp,1)/3e8;   % tau_cont; recompute to avoid closure issues
        nu  = 2*targets(pp,2)*77e9/3e8; % nu_hz; recompute (keeps consistent with top-level defaults)
        Hdelay_local = exp(-1j*2*pi * fvec * tau);
        tx_delay_local = ifft(TXF_local .* Hdelay_local, Lfft);
        tx_delay_local = tx_delay_local(1:Ltx_local);
        tx_delay_doppler_local = tx_delay_local .* exp(1j*2*pi*nu*tvec(1:Ltx_local));
        path_sig_local = hi * tx_delay_doppler_local;
        rx_nocp_sum_local = rx_nocp_sum_local + path_sig_local;
    end
    % attach CP and then receiver (remove CP and SFFT)
    rx_with_cp_local = [rx_nocp_sum_local(end-Mcp+1:end); rx_nocp_sum_local];
    rx_nocp_local = rx_with_cp_local(Mcp+1:end);
    Y_time_local = reshape(rx_nocp_local, N, M);
    Y_local = fft(Y_time_local, [], 1);
    y_vec = reshape(Y_local, [], 1);
end

% Build Heff by feeding unit impulses at each transmit DD-coefficient
% To ensure the same normalization and channel parameters are used inside the helper,
% we will call it with tx_norm_flag = false and do global normalizations externally if needed.
for tx_idx = 1:MN_tx
    e = zeros(MN_tx,1);
    e(tx_idx) = 1;
    y_e = txch_rx_map(e, N, M, Mcp, Lfft, fvec, tvec, targets, false);
    Heff(:, tx_idx) = y_e;
    if mod(tx_idx, 128) == 0
        fprintf('  Heff progress: %d / %d columns built\n', tx_idx, MN_tx);
    end
end
fprintf('Heff construction complete.\n');

% Estimate channel noise variance baseline for each SNR (we will add AWGN in evaluation)
% But Heff is built noiseless; for MMSE we need sigma2 (noise variance) corresponding to SNR of pilot.
% We'll define sigma2 from SNR value and pilot energy: pilot_energy = mean(|tx_cp|^2) (earlier)
pilot_energy = mean(abs(tx_cp).^2);
% For each SNR, perform nTrials Monte Carlo: add AWGN to y = Heff * x + w, apply MMSE, compute metrics.

BER_vs_SNR = zeros(1, nSNR);
RMSE_range_vs_SNR = zeros(1, nSNR);
RMSE_vel_vs_SNR = zeros(1, nSNR);

fprintf('\nStarting Monte Carlo sweep over SNRs...\n');
for iS = 1:nSNR
    SNRdB = SNRdB_list(iS);
    fprintf(' SNR = %d dB ... ', SNRdB);
    SNRlin = 10^(SNRdB/10);
    sigma2 = pilot_energy / SNRlin;
    nErr = 0;
    total_symbols = 0;
    rmse_range_accum = 0;
    rmse_vel_accum = 0;
    for trial = 1:nTrials
        % Generate random data grid Xgrid with same mask (data_mask) but new random QPSK
        xvec = zeros(MN,1);
        Xgrid_rand = zeros(N,M);
        rng(trial + iS*1000); % vary seed
        % Put pilot same as before
        Xgrid_rand(kp_m, lp_m) = Xp;
        % Fill data positions with random QPSK
        data_positions = find(data_mask);
        nDataPos = numel(data_positions);
        bits = randi([0 1], 2*nDataPos, 1);
        qpsk_sym = (1-2*bits(1:2:end)) + 1j*(1-2*bits(2:2:end));
        Xgrid_rand(data_mask) = qpsk_sym;
        xvec = reshape(Xgrid_rand, [], 1);
        
        % Create noiseless received vector y0 = Heff * xvec
        y0 = Heff * xvec;
        % Add AWGN
        w = sqrt(sigma2/2) * (randn(size(y0)) + 1j*randn(size(y0)));
        y = y0 + w;
        
        % MMSE equalizer: Xhat_vec = Heff' * inv(Heff*Heff' + sigma2 * I) * y
        % Compute inverse via Cholesky or backslash for numerical stability
        A = Heff * Heff' + sigma2 * eye(MN_rx);
        % Solve A * z = y  => z = A\y ; then Xhat = Heff' * z
        z = A \ y;
        Xhat_vec = Heff' * z;
        
        % Reshape to grid
        Xhat_grid = reshape(Xhat_vec, N, M);
        
        % Data symbol detection (QPSK): nearest quadrant
        recov_symbols = Xhat_grid(data_mask);
        orig_symbols = Xgrid_rand(data_mask);
        % Decision by sign of real/imag parts
        recov_dec = (sign(real(recov_symbols)) + 1)/-2 + 1j*(sign(imag(recov_symbols)) + 1)/-2; % not needed exactly; compute SER by signs
        % compute symbol errors
        errs = sum((sign(real(recov_symbols)) ~= sign(real(orig_symbols))) | (sign(imag(recov_symbols)) ~= sign(imag(orig_symbols))));
        nErr = nErr + errs;
        total_symbols = total_symbols + numel(orig_symbols);
        
        % Now do RVM-based range/velocity estimation on this noisy received vector y
        % Convert y to Y_grid form (inverse of reshape earlier)
        Y_hat_grid = reshape(y, N, M);
        % We performed earlier: RVM_complex = fft2(Y, Nfine, Mfine);
        RVMc = fft2(Y_hat_grid, Nfine, Mfine);
        RVM_loc = abs(fftshift(RVMc));
        % Find top Ltargets peaks
        RVM_copy_loc = RVM_loc;
        est_fine_indices_loc = zeros(Ltargets,2);
        for t = 1:Ltargets
            [~, idxp] = max(RVM_copy_loc(:));
            [r_idxp, c_idxp] = ind2sub(size(RVM_copy_loc), idxp);
            est_fine_indices_loc(t,:) = [c_idxp, r_idxp];
            RVM_copy_loc(r_idxp, c_idxp) = 0;
        end
        % Map fine indices to continuous bin coordinates (same mapping as earlier)
        li_hat_loc = zeros(Ltargets,1);
        ki_hat_loc = zeros(Ltargets,1);
        for t=1:Ltargets
            c_idxp = est_fine_indices_loc(t,1);
            r_idxp = est_fine_indices_loc(t,2);
            rel_colp = c_idxp - center_col;
            rel_rowp = r_idxp - center_row;
            li_coarse_p = mod(rel_colp/eps_theta, M);
            ki_coarse_p = mod(rel_rowp/eps_phi, N);
            li_hat_loc(t) = li_coarse_p;
            ki_hat_loc(t) = ki_coarse_p;
        end
        % Convert to physical range/velocity
        Ri_hat_loc = (li_hat_loc) * 3e8 ./ (2 * M * Delta_f);
        Vi_hat_loc = (ki_hat_loc) * 3e8 ./ (2 * fc * N * T);
        
        % Match estimated targets to ground truth by nearest-neighbor in (R,V) space
        true_R = targets(:,1);
        true_V = targets(:,2);
        est_R = Ri_hat_loc;
        est_V = Vi_hat_loc;
        matched_idx = zeros(Ltargets,1);
        squared_error_R = 0;
        squared_error_V = 0;
        % Greedy matching: for each true target, find closest estimated that isn't used
        used = false(Ltargets,1);
        for tt = 1:Ltargets
            dists = (est_R - true_R(tt)).^2 + (est_V - true_V(tt)).^2;
            dists(used) = Inf;
            [~, idxmin] = min(dists);
            used(idxmin) = true;
            matched_idx(tt) = idxmin;
            squared_error_R = squared_error_R + (est_R(idxmin) - true_R(tt))^2;
            squared_error_V = squared_error_V + (est_V(idxmin) - true_V(tt))^2;
        end
        rmse_range_accum = rmse_range_accum + sqrt(squared_error_R / Ltargets);
        rmse_vel_accum = rmse_vel_accum + sqrt(squared_error_V / Ltargets);
    end % trial loop
    
    BER_vs_SNR(iS) = nErr / total_symbols;
    RMSE_range_vs_SNR(iS) = rmse_range_accum / nTrials;
    RMSE_vel_vs_SNR(iS) = rmse_vel_accum / nTrials;
    fprintf(' done. BER=%.4e, RMSE_R=%.3f m, RMSE_V=%.3f m/s\n', BER_vs_SNR(iS), RMSE_range_vs_SNR(iS), RMSE_vel_vs_SNR(iS));
end % SNR loop

% Plot BER and RMSE curves
figure('Name','BER vs SNR');
semilogy(SNRdB_list, BER_vs_SNR, '-o', 'LineWidth', 1.8);
grid on; xlabel('SNR (dB)'); ylabel('BER'); title('BER vs SNR (MMSE Equalizer)');

figure('Name','RMSE vs SNR');
subplot(2,1,1);
plot(SNRdB_list, RMSE_range_vs_SNR, '-s', 'LineWidth', 1.6); grid on;
xlabel('SNR (dB)'); ylabel('RMSE Range (m)'); title('RMSE of Range vs SNR');
subplot(2,1,2);
plot(SNRdB_list, RMSE_vel_vs_SNR, '-d', 'LineWidth', 1.6); grid on;
xlabel('SNR (dB)'); ylabel('RMSE Velocity (m/s)'); title('RMSE of Velocity vs SNR');

fprintf('\nExtension finished. You now have BER and RMSE vs SNR curves.\n');