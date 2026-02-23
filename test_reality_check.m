% test_reality_check.m
% Verifies that Heff is non-trivial and tests performance under Noise.

clear; clc;

% 1. Load Results
if ~isfile('receiver_results.mat'), error('Run run_integration first.'); end
load('receiver_results.mat');
load('waveform_package.mat', 'X');

% 2. Check if Heff is just an Identity Matrix (Trivial)
diff_from_identity = norm(Heff - eye(size(Heff)), 'fro');
fprintf('--- DIAGNOSTIC 1: Is the Channel Real? ---\n');
if diff_from_identity < 1e-6
    fprintf(' WARNING: Heff is an Identity Matrix. The channel is doing nothing.\n');
else
    fprintf('PASS: Heff is complex (Diff = %.2f). The equalizer is working hard.\n', diff_from_identity);
end

% 3. Check Condition Number (Stability)
cond_num = cond(Heff);
fprintf('\n--- DIAGNOSTIC 2: Matrix Stability ---\n');
fprintf('Condition Number: %.2e\n', cond_num);
if cond_num > 1e5
    fprintf(' WARNING: Matrix is ill-conditioned. Noise will cause huge errors.\n');
else
    fprintf(' PASS: Matrix is stable.\n');
end

% 4. The "Noise Injection" Test
fprintf('\n--- DIAGNOSTIC 3: Adding Artificial Noise (20 dB SNR) ---\n');

% Generate Noise
sig_power = mean(abs(Y(:)).^2);
target_snr_db = 20; 
noise_power = sig_power / (10^(target_snr_db/10));
noise = sqrt(noise_power/2) * (randn(size(Y)) + 1j*randn(size(Y)));

% Add to Y
Y_noisy = Y + noise;

% Equalize Noisy Signal
Xhat_noisy_vec = mmse_equalize(Y_noisy, Heff, noise_power);
Xhat_noisy = reshape(Xhat_noisy_vec, size(X));

% Normalize (using Pilot Peak method)
[~, idxPeak] = max(abs(X(:)));
alpha = X(idxPeak) / Xhat_noisy(idxPeak);
Xhat_noisy = Xhat_noisy * alpha;

% Calculate Real SNR
maskData = abs(X) > 0; % Simplified mask
mse_noise = mean(abs(X(maskData) - Xhat_noisy(maskData)).^2);
snr_out = 10*log10( mean(abs(X(maskData)).^2) / mse_noise );

fprintf('Input SNR:  %.2f dB\n', target_snr_db);
fprintf('Output SNR: %.2f dB\n', snr_out);

if snr_out > 15 && snr_out < 25
    fprintf(' PASS: System behaves realistically under noise.\n');
elseif snr_out < 5
    fprintf(' FAIL: System collapsed under noise (Check MMSE regularization).\n');
else
    fprintf(' NOTE: Output SNR is unusually high/low.\n');
end

% 5. Visualize
figure('Name', 'Noise Test');
subplot(1,2,1); imagesc(abs(X)); title('Original X'); axis square;
subplot(1,2,2); imagesc(abs(Xhat_noisy)); title(['Recovered X (SNR ' num2str(snr_out, '%.1f') ' dB)']); axis square;