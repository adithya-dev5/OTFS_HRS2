function diagnostic_3()
% diagnostic_3.m
% Run a suite of automatic diagnostics on your OTFS chain and receiver_results.mat
% Fixed: Explicit variable initialization to prevent parser errors.

    close all; clc;
    
    % Initialize report container
    report = {};
    
    % --- Nested Helper: Push to Report ---
    function push(s)
        report{end+1} = s;
        fprintf('%s\n', s); % Also print to console
    end
    
    push('=== OTFS DIAGNOSTIC RUN ===');
    
    % Setup Paths (Manual addpath to be safe)
    if exist('src/receiver','dir'), addpath('src/receiver'); push('Added path: src/receiver'); end
    if exist('src/txr','dir'), addpath('src/txr'); push('Added path: src/txr'); end

    % 1) Try to load receiver_results.mat
    candidates = {'receiver_results.mat', fullfile(pwd,'receiver_results.mat'), '../receiver_results.mat'};
    s = struct(); loaded = false; loadedVars = {};
    
    for i=1:numel(candidates)
        f = candidates{i};
        if exist(f,'file')
            try
                s = load(f);
                loaded = true;
                loadedVars = fieldnames(s);
                push(['Loaded MAT: ' f]);
                push(['Variables: ' strjoin(loadedVars', ', ')]);
                break;
            catch ME
                push(['Failed reading ' f ': ' ME.message]);
            end
        end
    end
    
    if ~loaded
        push('No receiver_results.mat found. Aborting diagnostic.');
        return;
    end

    % --- 2) EXTRACT VARIABLES (Robustly) ---
    % Helper to get field or default
    function val = get_val(names, default_val)
        val = default_val;
        if ~iscell(names), names = {names}; end
        for k = 1:numel(names)
            if isfield(s, names{k})
                val = s.(names{k});
                push(['Mapped ' names{k} ' -> variable']);
                return;
            end
        end
    end

    % Initialize variables explicitly to avoid parser errors
    tx_sig = get_val({'tx_sig', 'tx'}, []);
    rx_int = get_val('rx_int', []);
    Y_int  = get_val({'Y_int', 'Y'}, []);
    X      = get_val('X', []);
    Xhat   = get_val('Xhat', []);
    Heff   = get_val('Heff', []);
    
    % Parameters
    M  = get_val('M', 32);
    N  = get_val('N', 32);
    lp = get_val('lp', 15);
    kp = get_val('kp', 15);
    
    push(sprintf('Using Parameters: M=%d, N=%d, lp=%d, kp=%d', M, N, lp, kp));

    % Masks
    [dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
    maskPilot_calc = (dd_L == lp) & (dd_K == kp);
    maskGuard_calc = ( abs(dd_L - lp) <= 2 ) | ( abs(dd_K - kp) <= 4 );
    maskData_calc = ~(maskGuard_calc | maskPilot_calc);

    maskPilot = get_val('maskPilot', maskPilot_calc);
    
    % Check Mask
    if isequal(maskPilot, maskPilot_calc)
        push('maskPilot matches calculated (OK).');
    else
        push('maskPilot in file differs from calculated defaults (Warning).');
    end

    % --- 3) RUN TESTS ---

    % Test 1: Mod/Demod Invertibility
    if exist('otfs_modulate','file') && exist('otfs_demodulate','file')
        push('Testing otfs_modulate -> otfs_demodulate invertibility...');
        X_test = zeros(N,M);
        X_test(kp+1, lp+1) = 1; % Pilot
        try
            tx_test = otfs_modulate(X_test, N, M);
            Y_test = otfs_demodulate(tx_test, 0, N, M);
            
            % Normalize for comparison
            if max(abs(Y_test(:))) > 0
                scale = max(abs(X_test(:))) / max(abs(Y_test(:)));
                Y_test = Y_test * scale;
            end
            
            reconstErr = norm(abs(Y_test(:)) - abs(X_test(:))) / max(1,norm(X_test(:)));
            push(sprintf('Reconstruction Error: %.3e', reconstErr));
            
            if reconstErr < 1e-6
                push('-> PASS: Invertibility OK.');
            else
                push('-> FAIL: Non-zero error. Check Transpose/Scaling.');
            end
        catch ME
            push(['Mod/Demod test failed: ' ME.message]);
        end
    else
        push('otfs_modulate/demodulate functions not found.');
    end

    % Test 2: Energy Check
    if ~isempty(X) && ~isempty(tx_sig)
        pX = mean(abs(X(:)).^2);
        pTx = mean(abs(tx_sig(:)).^2);
        push(sprintf('Mean Power: Grid=%.4e, TimeSig=%.4e', pX, pTx));
    end

    % Test 3: Pilot Location Check
    if ~isempty(Y_int)
        [~, idxLocal] = max(abs(Y_int(:)));
        [r_loc, c_loc] = ind2sub(size(Y_int), idxLocal);
        meas_k = r_loc - 1; meas_l = c_loc - 1;
        
        push(sprintf('Strongest Peak in Y: (Delay=%d, Doppler=%d)', meas_l, meas_k));
        push(sprintf('Expected Pilot:      (Delay=%d, Doppler=%d)', lp, kp));
        
        if meas_l == lp && meas_k == kp
            push('-> PASS: Pilot location matches.');
        else
            push('-> FAIL: Pilot location mismatch. Check Axes.');
        end
    else
        push('Y_int not found, skipping peak check.');
    end
    
    % Save Report
    fid = fopen('diagnostic_report.txt','w');
    if fid>0
        for ii=1:numel(report), fprintf(fid, '%s\n', report{ii}); end
        fclose(fid);
        fprintf('\nReport saved to diagnostic_report.txt\n');
    end
end