function RVM = compute_highres_RVM(Y, eps_phi, eps_theta)
% COMPUTE_HIGHRES_RVM: Interpolates the Delay-Doppler grid for visualization
% Y:         The coarse Delay-Doppler grid (from otfs_demodulate)
% eps_phi:   Oversampling factor for Delay (Rows)
% eps_theta: Oversampling factor for Doppler (Cols)

if nargin < 2, eps_phi = 4; end
if nargin < 3, eps_theta = 4; end

% Method 1: Sinc Interpolation (The DSP way)
% Transform to Dual Domain -> Zero Pad -> Inverse Transform
[M, N] = size(Y);
M_new = M * eps_phi;
N_new = N * eps_theta;

% 1. Transform DD -> Time-Frequency (Dual Domain)
% Note: Using fft2/ifft2 pair. Order matters less than consistency.
Y_dual = fft2(Y); 

% 2. Zero-Pad in the Dual Domain (High frequencies = Center of spectrum)
Y_dual_shifted = fftshift(Y_dual);
Y_dual_padded = zeros(M_new, N_new);

r_start = floor((M_new - M)/2) + 1;
c_start = floor((N_new - N)/2) + 1;
Y_dual_padded(r_start:r_start+M-1, c_start:c_start+N-1) = Y_dual_shifted;

% 3. Transform back to DD -> Interpolated
% Use ifft2 to reverse the earlier fft2
Y_interp = ifft2(ifftshift(Y_dual_padded));

% Normalize Magnitude
RVM = abs(Y_interp);
RVM = RVM / max(RVM(:) + eps);

% Optional: Use standard image interpolation (Simpler/Faster)
% [X, Z] = meshgrid(1:N, 1:M);
% [Xq, Zq] = meshgrid(linspace(1, N, N_new), linspace(1, M, M_new));
% RVM = interp2(X, Z, abs(Y), Xq, Zq, 'cubic');
% RVM = RVM / max(RVM(:));

end