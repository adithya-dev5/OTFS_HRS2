%% generate_annotated_rvm.m
% Generates a "Presentation Quality" RVM with Units and Arrows.
% Demonstrates exactly WHAT is being sensed (Range & Velocity).

clear; close all; clc;

% --- 1. SETUP PARAMETERS ---
if exist('src/txr','dir'), addpath('src/txr'); end
if exist('src/receiver','dir'), addpath('src/receiver'); end

% System Constants
M = 32; N = 32; 
Mcp = 32; 
fs = 15.36e6; 
fc = 77e9;         % 77 GHz (Automotive Radar Standard)
c_light = 3e8;

% Grid & Pilot
Xp = 10;
lp = 15; kp = 15;  % Pilot Center (0-based)
lTau = 2; kNu = 2;

% --- 2. CALCULATE PHYSICAL UNITS ---
dt = 1/fs;                  % Sample duration
df = fs/(M*N);              % Doppler bin resolution

% Radar Resolution Formulas
range_res = (c_light * dt) / 2;       % Meters per Sample
vel_res   = (c_light * df) / (2*fc);  % m/s per Bin

fprintf('System Resolution:\n');
fprintf('  Range:    %.2f meters/sample\n', range_res);
fprintf('  Velocity: %.2f m/s per bin\n', vel_res);

% --- 3. DEFINE TARGET (Physics) ---
% Target: 5.3 samples delay, 4.2 bins Doppler
target_delay_samples = 5.3;
target_doppler_bins  = 4.2;

t_frac = struct();
t_frac.gain = 0.8;
t_frac.tau = target_delay_samples * dt; 
t_frac.nu  = target_doppler_bins * df;

% Calculate Expected Physical Values
tgt_range_m = target_delay_samples * range_res;
tgt_vel_mps = target_doppler_bins * vel_res;

% --- 4. GENERATE SIGNAL & CHANNEL ---
% Build clean grid (Pilot + Data)
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskGuard = (abs(dd_L-lp)<=lTau) | (abs(dd_K-kp)<=2*kNu); 
maskPilot = (dd_L==lp) & (dd_K==kp);
maskData = ~(maskGuard | maskPilot);

rng(42); 
DataRef = qammod(randi([0 3], nnz(maskData), 1), 4, 'UnitAveragePower', true);
X = build_dd_grid(M, N, lp+1, kp+1, lTau, kNu, DataRef, Xp);

% Modulate & Channel
tx = otfs_modulate(X, N, M);
if Mcp>0, tx = [tx(end-Mcp+1:end); tx]; end

params = struct('M',M, 'N',N, 'Mcp',Mcp, 'fs',fs, 'Lfft',2048);
rx = apply_channel(tx, t_frac, params);

% Demodulate
Y = otfs_demodulate(rx, Mcp, N, M);

% --- 5. HIGH-RES INTERPOLATION ---
Q = 4; % 4x Interpolation Factor
Y_dual = fftshift(fft2(Y));
[Ny, Mx] = size(Y);
Y_pad = zeros(Q*Ny, Q*Mx);
r_s = floor((Q*Ny - Ny)/2) + 1;
c_s = floor((Q*Mx - Mx)/2) + 1;
Y_pad(r_s:r_s+Ny-1, c_s:c_s+Mx-1) = Y_dual;

RVM = abs(ifft2(ifftshift(Y_pad)));
RVM = RVM / max(RVM(:)); % Normalize

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
% --- 7. DRAW ARROWS & LABELS ---
hold on;

% 1. Mark Pilot (Origin)
plot(0, 0, 'wo', 'MarkerSize', 10, 'LineWidth', 2);
text(0, -5, 'Pilot (Ref)', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% 2. Mark Target
plot(tgt_range_m, tgt_vel_mps, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
text(tgt_range_m, tgt_vel_mps + 5, 'Target', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% 3. Draw Range Arrow (Horizontal)
quiver(0, 0, tgt_range_m, 0, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(tgt_range_m/2, -2, sprintf('%.1f m', tgt_range_m), 'Color', 'w', 'HorizontalAlignment', 'center');

% 4. Draw Velocity Arrow (Vertical, starting from Range point)
quiver(tgt_range_m, 0, 0, tgt_vel_mps, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text(tgt_range_m + 2, tgt_vel_mps/2, sprintf('%.1f m/s', tgt_vel_mps), 'Color', 'w');

% Zoom in on the action (Optional)
xlim([-10, tgt_range_m + 20]);
ylim([-20, tgt_vel_mps + 20]);

fprintf('Figure Generated.\n');