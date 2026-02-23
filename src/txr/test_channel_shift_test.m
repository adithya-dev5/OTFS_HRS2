%% test_channel_shift_test.m
clear; close all; clc;

% Check for receiver path
if exist('../receiver','dir'), addpath('../receiver'); end

% 1. Setup Config
% --- FIX 1: Set Mcp = 0 to avoid offset errors ---
% --- FIX 2: Set Xp = 10 to ensure Pilot is stronger than Data ---
cfg = struct('M',32,'N',32,'lp',15,'kp',15,'lTau',2,'kNu',2,'Xp',10+0j,'fs',15.36e6,'Mcp',0);

[X, tx, params, ~] = waveform_and_channel(cfg);
fs = params.fs;
N = params.N; 
M = params.M; 

% 2. Define Shifts
% --- DELAY: Choose 5 samples (Fractional Delay Test) ---
% In True OTFS (High-Res), 1 Sample Delay = 1 Delay Bin Shift.
d_samples = 5; 
tau = d_samples / fs;

% --- DOPPLER: Choose integer bin ---
k_bin = 4; 
nu = k_bin * (fs / (N*M)); 

fprintf('Testing: Delay = %d samples, Doppler Bin = %d\n', d_samples, k_bin);

% 3. Apply Channel
targets = struct('tau',tau,'nu',nu,'gain',1.0+0j);

if exist('apply_channel','file')
    rx = apply_channel(tx, targets, struct('fs',fs,'Lfft',params.Lfft,'t0',0));
else
    error('apply_channel.m not found');
end

% 4. Demodulate
% Pass 0 for Mcp (consistent with cfg.Mcp=0)
RX_dd = otfs_demodulate(rx, 0, N, M);

% 5. Find Peaks
[~, linidx] = max(abs(RX_dd(:)));
[rowIdx, colIdx] = ind2sub(size(RX_dd), linidx);

% Convert to 0-based
meas_doppler = rowIdx - 1;
meas_delay   = colIdx - 1;

% 6. Calculate Expected
lp = cfg.lp; kp = cfg.kp;

% EXPECTATION ALIGNMENT:
% Delay: 1 Bin = 1 Sample.
expected_delay = mod(lp + d_samples, M); 

% Doppler: Standard shift
expected_doppler = mod(kp + k_bin, N);

fprintf('------------------------------------------------\n');
fprintf('MEASURED: Delay (Col) = %d, Doppler (Row) = %d\n', meas_delay, meas_doppler);
fprintf('EXPECTED: Delay (Col) = %d, Doppler (Row) = %d\n', expected_delay, expected_doppler);
fprintf('------------------------------------------------\n');

% 7. Validation
delay_err = abs(meas_delay - expected_delay);
doppler_err = abs(meas_doppler - expected_doppler);

% Handle wrap-around diffs
if delay_err > M/2, delay_err = M - delay_err; end
if doppler_err > N/2, doppler_err = N - doppler_err; end

if delay_err <= 1 && doppler_err <= 1
    fprintf('PASS: Shifts match physics!\n');
 figure; imagesc(abs(RX_dd)); axis xy; colorbar; 
    title('RX DD Magnitude (Debug)');
    xlabel('Delay (Col)'); ylabel('Doppler (Row)');
    hold on; 
    plot(expected_delay+1, expected_doppler+1, 'ro', 'MarkerSize',10, 'LineWidth', 2);
    plot(meas_delay+1, meas_doppler+1, 'wx', 'MarkerSize',10, 'LineWidth', 2);
    legend('Expected','Measured');
else
    fprintf('FAIL: Mismatch larger than tolerance.\n');
    
    % Plot for visual debug
    figure; imagesc(abs(RX_dd)); axis xy; colorbar; 
    title('RX DD Magnitude (Debug)');
    xlabel('Delay (Col)'); ylabel('Doppler (Row)');
    hold on; 
    plot(expected_delay+1, expected_doppler+1, 'ro', 'MarkerSize',10, 'LineWidth', 2);
    plot(meas_delay+1, meas_doppler+1, 'wx', 'MarkerSize',10, 'LineWidth', 2);
    legend('Expected','Measured');
end