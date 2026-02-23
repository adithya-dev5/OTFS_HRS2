function rx = synth_channel_from_targets(tx, params)
% Simple narrowband channel: sum over targets of delayed & Doppler-shifted copies of tx.
% params.targets should be an array of structs having at least fields:
%   .li (integer delay index), .ki (integer doppler index), .gain (complex)
% or fields R (range index) and V (doppler index) etc.
% We assume sampling such that a shift by 'li' corresponds to integer sample delays.

tx = tx(:);
L = length(tx);
rx = zeros(size(tx));

t = (0:L-1).'; % sample index
if isfield(params,'targets')
    targ = params.targets;
    for i=1:numel(targ)
        % find integer delay & doppler fields robustly
        if isfield(targ(i),'li'), li = targ(i).li;
        elseif isfield(targ(i),'delay_idx'), li = targ(i).delay_idx;
        elseif isfield(targ(i),'R'), li = targ(i).R;
        else li = 0; end

        if isfield(targ(i),'ki'), ki = targ(i).ki;
        elseif isfield(targ(i),'doppler_idx'), ki = targ(i).doppler_idx;
        elseif isfield(targ(i),'V'), ki = targ(i).V;
        else ki = 0; end

        if isfield(targ(i),'gain'), g = targ(i).gain;
        elseif isfield(targ(i),'amplitude'), g = targ(i).amplitude;
        else g = 1; end

        % delay (circular for simplicity)
        shifted = circshift(tx, li);
        % Doppler modulation: complex exponential across samples
        doppler = exp(1j*2*pi*(ki/L)*t);
        rx = rx + g * shifted .* doppler;
    end
else
    error('No params.targets to synthesize channel.');
end

% add small noise if requested
if isfield(params,'noise_var')
    rx = rx + sqrt(params.noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));
end
end
