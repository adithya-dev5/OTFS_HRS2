function peaks = find_rvm_peaks(RVM, Ltargets, cfg)
% FIND_RVM_PEAKS Canonical peak finder with explicit row/col semantics.
assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
if nargin < 2, Ltargets = 1; end
R = RVM;
peaks = zeros(Ltargets,2); % [doppler_row, delay_col]
mask = true(size(R));
for i=1:Ltargets
    Rmasked = R .* mask;
    [~, idx] = max(Rmasked(:));
    [r,c] = ind2sub(size(Rmasked), idx);
    peaks(i,:) = [r,c];
    rad = 6;
    rr = max(1, r-rad):min(size(R,1), r+rad);
    cc = max(1, c-rad):min(size(R,2), c+rad);
    mask(rr,cc) = false;
end
end
