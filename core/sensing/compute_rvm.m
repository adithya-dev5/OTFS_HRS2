function RVM = compute_rvm(Y, cfg, q_doppler, q_delay)
% COMPUTE_RVM Canonical RVM interpolation under explicit DD axis convention.
assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
if nargin < 3, q_doppler = 4; end
if nargin < 4, q_delay = 4; end

Y_dual = fftshift(fft2(Y));
[Ny, Mx] = size(Y);
Y_padded = zeros(q_doppler*Ny, q_delay*Mx);
r_start = floor((q_doppler*Ny - Ny)/2) + 1;
c_start = floor((q_delay*Mx - Mx)/2) + 1;
Y_padded(r_start:r_start+Ny-1, c_start:c_start+Mx-1) = Y_dual;
RVM = abs(ifft2(ifftshift(Y_padded)));
RVM = RVM / max(RVM(:) + eps);
end
