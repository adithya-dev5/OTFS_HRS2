function tx = otfs_modulate(X, N, M)
% OTFS modulate: High-Resolution (Delay=FastTime)
% LINEAR VERSION: No dynamic normalization.
% Input X: Rows=Doppler (N), Cols=Delay (M)

[Na, Ma] = size(X);
if Na ~= N || Ma ~= M
   if Na == M && Ma == N, X = X.'; else, error('X size mismatch. Expected %dx%d', N, M); end
end

% --- STEP 1: ISFFT Partial (Doppler -> Coarse Time) ---
% IFFT along Doppler (Rows)
X_time_delay = ifft(X, [], 1);

% --- STEP 2: Transpose (Align Delay with Fast Time) ---
X_fast_slow = X_time_delay.';

% --- STEP 3: Serialize ---
tx = X_fast_slow(:);

% --- NO NORMALIZATION (Crucial for Heff Linearity) ---
% tx = tx / std(tx);  <-- THIS LINE CAUSED THE BER BUG
end