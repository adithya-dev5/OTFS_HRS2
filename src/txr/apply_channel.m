function rx = apply_channel(tx, targets, params)
% APPLY_CHANNEL: High-Fidelity Physics (Fractional Delay + Doppler)
% FIXED: Phase wrapping for fractional delays.

tx = tx(:);
if ~isfield(params,'fs'), error('params.fs is required'); end
fs = params.fs;

% Determine FFT size (must cover full signal)
if ~isfield(params,'Lfft') || isempty(params.Lfft)
    Lfft = 2^nextpow2(max(4*length(tx), length(tx)));
else
    Lfft = max(params.Lfft, length(tx));
end

t0 = 0; if isfield(params,'t0'), t0 = params.t0; end
tvec = t0 + (0:length(tx)-1).' / fs;

% --- FIX: CORRECT FREQUENCY VECTOR ---
% Map upper half [N/2 ... N-1] to negative frequencies [-N/2 ... -1]
% This ensures exp(-j*w*tau) is continuous at Nyquist.
df = fs / Lfft;
k_idx = (0:Lfft-1).';
k_idx(k_idx >= Lfft/2) = k_idx(k_idx >= Lfft/2) - Lfft;
fvec = k_idx * df;
% -------------------------------------

TX_F = fft(tx, Lfft);
rx = zeros(size(tx));

for p = 1:numel(targets)
    tg = targets(p);
    
    % Extract parameters robustly
    if isfield(tg,'tau'), tau = tg.tau; elseif isfield(tg,'delay'), tau = tg.delay; else, tau=0; end
    if isfield(tg,'nu'), nu = tg.nu; elseif isfield(tg,'doppler'), nu = tg.doppler; else, nu=0; end
    if isfield(tg,'gain'), gain = tg.gain; else, gain=1; end
    
    % 1. Apply Delay (Phase Ramp in Freq Domain)
    phase_delay = exp(-1j*2*pi * fvec * tau);
    X_shifted = TX_F .* phase_delay;
    
    % 2. Transform back to Time
    x_delayed_full = ifft(X_shifted, Lfft);
    x_delayed = x_delayed_full(1:length(tx));
    
    % 3. Apply Doppler (Phase Ramp in Time Domain)
    doppler_phasor = exp(1j*2*pi * nu * tvec);
    x_doppler = x_delayed .* doppler_phasor;
    
    rx = rx + gain * x_doppler;
end
end