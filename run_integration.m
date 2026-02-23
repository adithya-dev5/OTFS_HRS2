% run_integration.m (Final Physics-Compliant Version)
% Trusted Physics: Uses Heff to correct the channel, no manual shifting.

clear; close all; clc;

% Ensure paths
if exist('src/receiver','dir'), addpath('src/receiver'); end
if exist('src/txr','dir'), addpath('src/txr'); end

% -------------------------
% 1. Load waveform package
% -------------------------
pkgFilename = 'waveform_package.mat';
if ~isfile(pkgFilename)
    error('Cannot find %s. Run regenerate_waveform.m first.', pkgFilename);
end
S = load(pkgFilename);

% Extract variables
if isfield(S,'tx_signal'), tx = S.tx_signal(:); elseif isfield(S,'tx'), tx = S.tx(:); else tx = []; end
if isfield(S,'X'), X = S.X; else X = []; end
if isfield(S,'params'), params = S.params; else params = struct(); end

% PATCH: Targets
if ~isfield(params, 'targets')
    if isfield(S, 'targets'), params.targets = S.targets;
    else, params.targets = struct('tau',0,'nu',0,'gain',0); end
end
% Inference
if ~isfield(params,'M') || ~isfield(params,'N')
    if ~isempty(X), [params.N, params.M] = size(X);
    else, params.M = 32; params.N = 32; end
end
if ~isfield(params,'Mcp'), params.Mcp = 0; end
if ~isfield(params,'Lfft'), params.Lfft = max(1024, numel(tx)); end

% -------------------------
% 2. Generate Rx (Physics)
% -------------------------
fprintf('Generating Channel Response...\n');
if exist('apply_channel','file')
    rx = apply_channel(tx, params.targets, params);
else
    rx = synth_channel_from_targets(tx, params);
end

% Add Noise (30 dB SNR)
rx_power = mean(abs(rx).^2);
noise_var = rx_power / (10^(30/10)); 
rx = rx + sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));

% Add CP handling
if params.Mcp > 0, rx = [rx(end-params.Mcp+1:end); rx]; end

% -------------------------
% 3. Demodulate
% -------------------------
Y = otfs_demodulate(rx, params.Mcp, params.N, params.M);

% -------------------------
% 4. Equalization (MMSE)
% -------------------------
fprintf('Building Heff (This captures the physics)...\n');
if exist('build_Heff','file')
    Heff = build_Heff(params.N, params.M, params, tx);
else
    Heff = eye(params.N * params.M);
    warning('build_Heff not found! Equalizer will fail.');
end

fprintf('Running MMSE Equalizer...\n');
Xhat_vec = mmse_equalize(Y, Heff, noise_var);
Xhat = reshape(Xhat_vec, params.N, params.M);

% Normalize Scale (Pilot-based)
if ~isempty(X)
    [~, idx] = max(abs(X(:)));
    scale = X(idx) / Xhat(idx);
    Xhat = Xhat * scale;
end

% -------------------------
% 5. Visualization & Results
% -------------------------

% --- STEP 5a: Calculate Metrics FIRST ---
snr_out = 0; % Default
if ~isempty(X)
    mask = abs(X) > 0.1; % Exclude pilot/guard for SNR stats
    if any(mask(:))
        mse = mean(abs(X(mask) - Xhat(mask)).^2);
        snr_out = 10*log10(mean(abs(X(mask)).^2) / mse);
        
        fprintf('\n======================================\n');
        fprintf('FINAL PHYSICS RESULTS:\n');
        fprintf('Recovered SNR: %.2f dB (Should be > 20 dB)\n', snr_out);
        fprintf('======================================\n');
    end
end

% --- STEP 5b: Plotting ---

% 1. Sensing RVM (High-Res Interpolation)
if exist('compute_highres_RVM', 'file')
    RVM = compute_highres_RVM(Y, 4, 4);
else
    % Manual Interpolation (FFT -> Pad -> IFFT)
    Y_dual = fftshift(fft2(Y));
    [Ny, Mx] = size(Y);
    Y_padded = zeros(4*Ny, 4*Mx);
    r_start = floor((4*Ny - Ny)/2) + 1;
    c_start = floor((4*Mx - Mx)/2) + 1;
    Y_padded(r_start:r_start+Ny-1, c_start:c_start+Mx-1) = Y_dual;
    RVM = abs(ifft2(ifftshift(Y_padded)));
end
RVM = RVM / max(RVM(:));

f1 = figure('Name', 'Sensing RVM', 'Color', 'w');
doppler_axis = linspace(0, params.N-1, size(RVM,1));
delay_axis   = linspace(0, params.M-1, size(RVM,2));
imagesc(delay_axis, doppler_axis, RVM); 
axis xy; colormap(jet); colorbar;
title('Sensing: Range-Velocity Map');
xlabel('Delay Bin (Fast Time)'); ylabel('Doppler Bin (Slow Time)');

% Highlight Target
hold on;
plot(15, 15, 'wo', 'MarkerSize', 10, 'LineWidth', 2); % Pilot
if isfield(params, 'lp')
    % Expected location (+5, +8)
    t_del = mod(params.lp + 5, params.M);
    t_dop = mod(params.kp + 8, params.N);
    plot(t_del, t_dop, 'wx', 'MarkerSize', 15, 'LineWidth', 2);
end
drawnow;

% 2. Comms Constellation
f2 = figure('Name', 'Comms Constellation', 'Color', 'w');
subplot(1,2,1); 
imagesc(abs(X)); title('TX Grid Magnitude'); axis xy; colorbar;
xlabel('Delay'); ylabel('Doppler');

subplot(1,2,2); 
if ~isempty(X)
    maskData = (abs(X) > 0.1) & (abs(X) < 9); % Data only
    plot(X(maskData), 'bo', 'MarkerSize', 4); hold on;
    plot(Xhat(maskData), 'rx', 'MarkerSize', 4); 
    legend('Transmitted', 'Recovered');
    title(sprintf('Constellation (SNR: %.2f dB)', snr_out)); 
    axis equal; grid on; box on;
end
drawnow;

save('receiver_results.mat', 'Y', 'RVM', 'Heff', 'Xhat'); % Force display

% Metrics
if ~isempty(X)
    mask = abs(X) > 0;
    mse = mean(abs(X(mask) - Xhat(mask)).^2);
    snr_out = 10*log10(mean(abs(X(mask)).^2) / mse);
    
    fprintf('\n======================================\n');
    fprintf('FINAL PHYSICS RESULTS:\n');
    fprintf('Recovered SNR: %.2f dB (Should be > 20 dB)\n', snr_out);
    fprintf('======================================\n');
else
    fprintf('No Ground Truth X found. Skipping SNR.\n');
end

save('receiver_results.mat', 'Y', 'RVM', 'Heff', 'Xhat');