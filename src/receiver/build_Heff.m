function Heff = build_Heff(N, M, params, ~)
% BUILD_HEFF Constructs the effective channel matrix Heff.
% FIXED: Loop order matches MATLAB's column-major X(:) standard.

if ~exist('otfs_modulate','file') || ~exist('otfs_demodulate','file')
    error('Modulator/Demodulator not found on path.');
end

MN = M * N;
Heff = zeros(MN, MN); 

Mcp = 0;
if isfield(params, 'Mcp'), Mcp = params.Mcp; end

% --- DIAGNOSTIC PRINT ---
% We print this to CONFIRM the new file is active
fprintf('Building Heff (Order: Delay-Outer, Doppler-Inner)...');

col_idx = 0;

% --- FIX: Match MATLAB's X(:) Vectorization Order ---
% X(:) stacks columns (Delay bins). 
% So we must iterate through Columns (m) slowly, and Rows (n) quickly.

for m = 1:M        % Delay Index (Cols of X) - OUTER LOOP
    for n = 1:N    % Doppler Index (Rows of X) - INNER LOOP
        col_idx = col_idx + 1;
        
        % 1. Create Impulse at (n,m)
        X_imp = zeros(N, M);
        X_imp(n, m) = 1;
        
        % 2. Modulate
        tx_imp = otfs_modulate(X_imp, N, M);
        if Mcp > 0, tx_imp = [tx_imp(end-Mcp+1:end); tx_imp]; end
        
        % 3. Channel
        if exist('apply_channel','file')
            rx_imp = apply_channel(tx_imp, params.targets, params);
        else
            rx_imp = synth_channel_from_targets(tx_imp, params);
        end
        
        % 4. Demodulate
        Y_imp = otfs_demodulate(rx_imp, Mcp, N, M);
        
        % 5. Store Vectorized Response
        Heff(:, col_idx) = Y_imp(:);
    end
    % Progress dot every few columns
    if mod(col_idx, 50) == 0, fprintf('.'); end
end
fprintf(' Done.\n');

end