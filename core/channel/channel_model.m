function rx = channel_model(tx, cfg, mode)
% CHANNEL_MODEL Unified channel interface.
% mode: "synthetic_integer" | "synthetic_fractional" | "backscatter"
if nargin < 3 || strlength(mode) == 0
    mode = cfg.channel.mode;
end

switch string(mode)
    case "synthetic_integer"
        rx = synthetic_integer_channel(tx, cfg);
    case "synthetic_fractional"
        rx = synthetic_fractional_channel(tx, cfg);
    case "backscatter"
        error('backscatter mode is a placeholder stub (not implemented).');
    otherwise
        error('Unknown channel mode: %s', mode);
end
end

function rx = synthetic_integer_channel(tx, cfg)
% Preserves prior integer-bin synthetic channel behavior.
tx = tx(:); L = length(tx); t = (0:L-1).';
rx = zeros(size(tx));
targ = cfg.targets;
for i = 1:numel(targ)
    if isfield(targ(i),'li'), li = targ(i).li;
    elseif isfield(targ(i),'delay_idx'), li = targ(i).delay_idx;
    elseif isfield(targ(i),'R'), li = targ(i).R;
    else, li = round(targ(i).tau * cfg.fs); end

    if isfield(targ(i),'ki'), ki = targ(i).ki;
    elseif isfield(targ(i),'doppler_idx'), ki = targ(i).doppler_idx;
    elseif isfield(targ(i),'V'), ki = targ(i).V;
    else, ki = round(targ(i).nu / (cfg.fs/(cfg.M*cfg.N))); end

    if isfield(targ(i),'gain'), g = targ(i).gain; else, g = 1; end

    shifted = circshift(tx, li);
    doppler = exp(1j*2*pi*(ki/L)*t);
    rx = rx + g * shifted .* doppler;
end
end

function rx = synthetic_fractional_channel(tx, cfg)
% Preserves prior fractional delay + continuous doppler channel behavior.
tx = tx(:);
fs = cfg.fs;
if isfield(cfg,'channel') && isfield(cfg.channel,'Lfft') && ~isempty(cfg.channel.Lfft)
    Lfft = max(cfg.channel.Lfft, length(tx));
else
    Lfft = 2^nextpow2(max(4*length(tx), length(tx)));
end

t0 = 0;
if isfield(cfg,'t0'), t0 = cfg.t0; end
tvec = t0 + (0:length(tx)-1).' / fs;

df = fs / Lfft;
k_idx = (0:Lfft-1).';
k_idx(k_idx >= Lfft/2) = k_idx(k_idx >= Lfft/2) - Lfft;
fvec = k_idx * df;

TX_F = fft(tx, Lfft);
rx = zeros(size(tx));
for p = 1:numel(cfg.targets)
    tg = cfg.targets(p);
    if isfield(tg,'tau'), tau = tg.tau; else, tau = 0; end
    if isfield(tg,'nu'), nu = tg.nu; else, nu = 0; end
    if isfield(tg,'gain'), gain = tg.gain; else, gain = 1; end

    phase_delay = exp(-1j*2*pi * fvec * tau);
    X_shifted = TX_F .* phase_delay;
    x_delayed = ifft(X_shifted, Lfft);
    x_delayed = x_delayed(1:length(tx));

    doppler_phasor = exp(1j*2*pi * nu * tvec);
    rx = rx + gain * (x_delayed .* doppler_phasor);
end
end
