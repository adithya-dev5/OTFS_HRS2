function Y = otfs_demodulate(rx_cp, Mcp, N, M)
    % 1. Remove CP and Vectorize
    rx_cp = rx_cp(:);
    if Mcp > 0, rx_noCP = rx_cp(Mcp+1:end); else, rx_noCP = rx_cp; end
    
    % Safety Truncation/Padding
    len = N * M;
    if numel(rx_noCP) < len, rx_noCP = [rx_noCP; zeros(len - numel(rx_noCP),1)];
    elseif numel(rx_noCP) > len, rx_noCP = rx_noCP(1:len); end
    
    % --- STEP 1: Reshape to Fast/Slow Time ---
    % We modulated M samples per block (Fast Time = Delay).
    % So we reshape to [M, N].
    X_fast_slow = reshape(rx_noCP, [M, N]);
    
    % --- STEP 2: Transpose back ---
    % Result: Rows=CoarseTime, Cols=Delay
    X_time_delay = X_fast_slow.';
    
    % --- STEP 3: SFFT Partial (Coarse Time -> Doppler) ---
    % FFT along Coarse Time (Rows)
    % Result: Rows=Doppler, Cols=Delay
    Y = fft(X_time_delay, [], 1);
end