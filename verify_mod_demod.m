% verify_mod_demod_v2.m
clear; clc; close all;

% 1. Setup Parameters
N = 32; % Doppler bins
M = 32; % Delay bins
Mcp = 0; 

% 2. Create an ASYMMETRIC Pattern (Letter 'L')
X_tx = zeros(N, M);
X_tx(5:10, 5) = 1;      % Vertical line 
X_tx(10, 5:10) = 1;     % Horizontal line

% 3. Modulate
tx_sig = otfs_modulate(X_tx, N, M);

% 4. Demodulate (Perfect Channel)
rx_sig = tx_sig; 
Y_rx_raw = otfs_demodulate(rx_sig, Mcp, N, M);

% --- CRITICAL FIX: Normalize Y to match X's Scale ---
% Because otfs_modulate boosts power, we must scale Y back down to compare.
scale_factor = max(abs(Y_rx_raw(:))) / max(abs(X_tx(:)));
Y_rx = Y_rx_raw / scale_factor;

fprintf('Detected Scale Factor: %.4f (Compensated)\n', scale_factor);

% 5. Visualization
figure('Position',[100 100 1000 400]);

subplot(1,3,1);
imagesc(abs(X_tx)); axis xy; title('TX (Original)');
xlabel('Delay (M)'); ylabel('Doppler (N)');
colorbar;

subplot(1,3,2);
imagesc(abs(Y_rx)); axis xy; title('RX (Scaled)');
xlabel('Delay (M)'); ylabel('Doppler (N)');
colorbar;

subplot(1,3,3);
% Calculate difference
diff_grid = abs(X_tx - Y_rx);
imagesc(diff_grid); axis xy; title('Difference');
colorbar;

% 6. Robust Diagnosis
% Flatten everything to vectors for safe comparison
vec_X = X_tx(:);
vec_Y = Y_rx(:);
vec_Y_trans = reshape(Y_rx.', [], 1); % Vector of Transposed Y

err_direct = norm(abs(vec_X - vec_Y));
err_trans  = norm(abs(vec_X - vec_Y_trans));

fprintf('--------------------------------------------------\n');
fprintf('Direct Error:     %.4f\n', err_direct);
fprintf('Transpose Error:  %.4f\n', err_trans);

if err_direct < 1e-1
    fprintf('RESULT: MATCH! (Mod/Demod are aligned)\n');
elseif err_trans < 1e-1
    fprintf("RESULT: ️ TRANSPOSED!\nFIX: Add .\' to the end of otfs_demodulate.m\n");
else
    fprintf('RESULT: MISMATCH. \nCheck if otfs_demodulate uses "fft(..., 2)". It should NOT.\n');
end
fprintf('--------------------------------------------------\n');