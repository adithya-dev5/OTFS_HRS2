function metrics = evaluate_ISAC(X, Xhat, params, peaks)
% Robust evaluation routine for ISAC pipeline.
% - X, Xhat: MN grids (can be empty)
% - params: struct loaded from waveform_package.mat (may contain params.targets)
% - peaks: Lx2 peak indices returned by find_RVM_peaks (row,col)
%
% Returns struct with:
%  .SER
%  .peaks
%  .RMSE_range_idx
%  .RMSE_doppler_idx
%  .notes

metrics = struct();
metrics.notes = {};

%% 1) Symbol-error-rate (SER) over data bins (if X present)
try
    if ~isempty(X)
        % Use data mask if available, otherwise non-zero X entries
        if isfield(params,'maskData') && ~isempty(params.maskData)
            data_idx = params.maskData;
        else
            data_idx = (abs(X) > 0);
        end

        if any(data_idx(:))
            txsym = X(data_idx);
            rxsym = Xhat(data_idx);

            % Simple decision rule: nearest quadrant (works for BPSK/QPSK-ish)
            tx_dec = sign(real(txsym)) + 1j*sign(imag(txsym));
            rx_dec = sign(real(rxsym)) + 1j*sign(imag(rxsym));
            ser = sum(tx_dec ~= rx_dec) / numel(tx_dec);
        else
            ser = NaN;
            metrics.notes{end+1} = 'No data symbols found for SER calculation.';
        end
    else
        ser = NaN;
        metrics.notes{end+1} = 'X (reference DD grid) is empty; SER not computed.';
    end
catch ME
    ser = NaN;
    metrics.notes{end+1} = ['SER computation failed: ' ME.message];
end
metrics.SER = ser;

%% 2) Peak bookkeeping
metrics.peaks = peaks;

%% 3) Try to extract true DD indices from params.targets robustly
true_li = [];
true_ki = [];

if isfield(params,'targets') && ~isempty(params.targets)
    targ = params.targets;
    Ltrue = numel(targ);
    true_li = nan(Ltrue,1);
    true_ki = nan(Ltrue,1);

    % helper flags for conversion
    hasDeltaF = isfield(params,'Delta_f') && ~isempty(params.Delta_f);
    hasM = isfield(params,'M') && ~isempty(params.M);
    hasN = isfield(params,'N') && ~isempty(params.N);
    hasT = isfield(params,'T') && ~isempty(params.T);
    hasfc = isfield(params,'fc') || isfield(params,'f_c') || isfield(params,'carrier_freq');

    % canonical parameter names
    if ~hasT && hasDeltaF
        Tval = 1/params.Delta_f;
    elseif hasT
        Tval = params.T;
    else
        Tval = [];
    end

    if hasfc
        if isfield(params,'fc'), fc = params.fc;
        elseif isfield(params,'f_c'), fc = params.f_c;
        else fc = params.carrier_freq;
        end
    else
        fc = [];
    end

    c0 = 3e8; % speed of light

    for ii=1:Ltrue
        s = targ(ii);

        % 1) delay index candidates
        if isfield(s,'li'), true_li(ii) = s.li; end
        if isnan(true_li(ii)) && isfield(s,'delay_idx'), true_li(ii) = s.delay_idx; end
        if isnan(true_li(ii)) && isfield(s,'delayIndex'), true_li(ii) = s.delayIndex; end
        if isnan(true_li(ii)) && isfield(s,'delay'), true_li(ii) = s.delay; end

        % 2) doppler index candidates
        if isfield(s,'ki'), true_ki(ii) = s.ki; end
        if isnan(true_ki(ii)) && isfield(s,'doppler_idx'), true_ki(ii) = s.doppler_idx; end
        if isnan(true_ki(ii)) && isfield(s,'dopplerIndex'), true_ki(ii) = s.dopplerIndex; end
        if isnan(true_ki(ii)) && isfield(s,'doppler'), true_ki(ii) = s.doppler; end

        % 3) try conversion from physical R (meters) -> li, using:
        %    tau = 2*R/c0 ; li ~= tau * M * Delta_f  (based on li/(M*Delta_f) ~= round-trip delay)
        if isnan(true_li(ii))
            % look for range fields
            if isfield(s,'li_m') && ~isempty(s.li_m)
                true_li(ii) = s.li_m; % already an index-like value
            elseif isfield(s,'R') && ~isempty(s.R) && hasDeltaF && hasM
                Rm = s.R;
                tau = 2*Rm / c0;
                true_li(ii) = round(tau * params.M * params.Delta_f);
                metrics.notes{end+1} = sprintf('Converted R->li for target %d using params.M and params.Delta_f',ii);
            elseif isfield(s,'range_m') && ~isempty(s.range_m) && hasDeltaF && hasM
                Rm = s.range_m;
                tau = 2*Rm / c0;
                true_li(ii) = round(tau * params.M * params.Delta_f);
                metrics.notes{end+1} = sprintf('Converted range_m->li for target %d',ii);
            end
        end

        % 4) try conversion from velocity V (m/s) -> ki, using:
        %    nu = 2*V*fc / c0 ; ki ~= nu * N * T
        if isnan(true_ki(ii))
            if isfield(s,'V') && ~isempty(s.V) && ~isempty(fc) && hasN && ~isempty(Tval)
                Vi = s.V;
                nu = 2 * Vi * fc / c0;
                true_ki(ii) = round(nu * params.N * Tval);
                metrics.notes{end+1} = sprintf('Converted V->ki for target %d using params.N, T, fc',ii);
            elseif isfield(s,'vel_mps') && ~isempty(s.vel_mps) && ~isempty(fc) && hasN && ~isempty(Tval)
                Vi = s.vel_mps;
                nu = 2 * Vi * fc / c0;
                true_ki(ii) = round(nu * params.N * Tval);
                metrics.notes{end+1} = sprintf('Converted vel_mps->ki for target %d',ii);
            end
        end

        % If some fields are strings, attempt numeric conversion
        if ~isnumeric(true_li(ii)) && ~isempty(true_li(ii))
            try true_li(ii) = str2double(true_li(ii)); catch, true_li(ii)=NaN; end
        end
        if ~isnumeric(true_ki(ii)) && ~isempty(true_ki(ii))
            try true_ki(ii) = str2double(true_ki(ii)); catch, true_ki(ii)=NaN; end
        end
    end
else
    metrics.notes{end+1} = 'params.targets not present; skipping RMSE target evaluation.';
end
%% 4) Compute RMSE between found peaks and true indices
try
    if ~isempty(true_li) && ~all(isnan(true_li)) && ~isempty(peaks)
        % peaks is [row, col]
        % Since grid is (N x M) = (Doppler x Delay):
        % Row (dim 1) = Doppler (ki)
        % Col (dim 2) = Delay (li)
        
        Lp = size(peaks,1);
        Ltrue = numel(true_li);
        Lmin = min([Lp, Ltrue]);
        
        if Lmin >= 1
            % --- CORRECTION HERE ---
            est_ki = peaks(1:Lmin, 1); % Row is Doppler
            est_li = peaks(1:Lmin, 2); % Col is Delay
            % -----------------------
            
            valid_li = ~isnan(true_li(1:Lmin));
            valid_ki = ~isnan(true_ki(1:Lmin));
            
            if any(valid_li)
                metrics.RMSE_range_idx = sqrt(mean((est_li(valid_li) - true_li(valid_li)).^2));
            else
                metrics.RMSE_range_idx = NaN;
                metrics.notes{end+1} = 'No valid true li values.';
            end
            
            if any(valid_ki)
                metrics.RMSE_doppler_idx = sqrt(mean((est_ki(valid_ki) - true_ki(valid_ki)).^2));
            else
                metrics.RMSE_doppler_idx = NaN;
                metrics.notes{end+1} = 'No valid true ki values.';
            end
        else
            % ... [Keep existing else block] ...
            metrics.RMSE_range_idx = NaN;
            metrics.RMSE_doppler_idx = NaN;
        end
    else
        % ... [Keep existing else block] ...
        metrics.RMSE_range_idx = NaN;
        metrics.RMSE_doppler_idx = NaN;
    end
catch ME
    metrics.RMSE_range_idx = NaN;
    metrics.RMSE_doppler_idx = NaN;
    metrics.notes{end+1} = ['RMSE computation failed: ' ME.message];
end

%% 4) Compute RMSE between found peaks and true indices (if possible)
try
    if ~isempty(true_li) && ~all(isnan(true_li)) && ~isempty(peaks)
        % peaks is Lx2: [row, col] likely corresponds to [li, ki] mapping depending on RVM orientation.
        % Here we assume peaks(:,1) -> row -> li index and peaks(:,2) -> col -> ki index
        Lp = size(peaks,1);
        Ltrue = numel(true_li);
        Lmin = min([Lp, Ltrue]);
        if Lmin >= 1
            est_li = peaks(1:Lmin,1);
            est_ki = peaks(1:Lmin,2);
            valid_li = ~isnan(true_li(1:Lmin));
            valid_ki = ~isnan(true_ki(1:Lmin));
            if any(valid_li)
                metrics.RMSE_range_idx = sqrt(mean((est_li(valid_li) - true_li(valid_li)).^2));
            else
                metrics.RMSE_range_idx = NaN;
                metrics.notes{end+1} = 'No valid true li values for RMSE computation.';
            end
            if any(valid_ki)
                metrics.RMSE_doppler_idx = sqrt(mean((est_ki(valid_ki) - true_ki(valid_ki)).^2));
            else
                metrics.RMSE_doppler_idx = NaN;
                metrics.notes{end+1} = 'No valid true ki values for RMSE computation.';
            end
        else
            metrics.RMSE_range_idx = NaN;
            metrics.RMSE_doppler_idx = NaN;
            metrics.notes{end+1} = 'Not enough peaks or targets to compute RMSE.';
        end
    else
        metrics.RMSE_range_idx = NaN;
        metrics.RMSE_doppler_idx = NaN;
        metrics.notes{end+1} = 'true_li/true_ki empty or peaks empty; RMSE not computed.';
    end
catch ME
    metrics.RMSE_range_idx = NaN;
    metrics.RMSE_doppler_idx = NaN;
    metrics.notes{end+1} = ['RMSE computation failed: ' ME.message];
end

end
