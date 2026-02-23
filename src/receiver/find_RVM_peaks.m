function peaks = find_RVM_peaks(RVM, Ltargets)
% picks top Ltargets peaks with simple non-max suppression
if nargin < 2, Ltargets = 1; end

R = RVM;
peaks = zeros(Ltargets,2);

% Simple iterative greedy peak picking with local suppression
mask = true(size(R));
for i=1:Ltargets
    Rmasked = R .* mask;
    [~, idx] = max(Rmasked(:));
    [r,c] = ind2sub(size(Rmasked), idx);
    peaks(i,:) = [r,c];
    % suppress neighborhood
    rad = 6; % suppression radius (pixels) - tune if needed
    rr = max(1, r-rad):min(size(R,1), r+rad);
    cc = max(1, c-rad):min(size(R,2), c+rad);
    mask(rr,cc) = false;
end
end
