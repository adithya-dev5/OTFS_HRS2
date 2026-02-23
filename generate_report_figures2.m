function generate_final_report_figures()
%% generate_final_report_figures.m
% Generates Figures 5.1 to 5.7 for the Final Report.
% FIXED: Self-contained Grid Logic to prevent index mismatch.

    clearvars -except generate_final_report_figures; 
    close all; clc;

    % --- 0. PATHS & SETUP ---
    if exist('src/txr','dir'), addpath('src/txr'); end
    if exist('src/receiver','dir'), addpath('src/receiver'); end

    % System Constants
    M = 32; N = 32; 
    Mcp = 32;          % CP Enabled
    fs = 15.36e6; 
    Xp = 10;           % Strong Pilot
    lp = 15; kp = 15;  % Center (0-based)
    lTau = 2; kNu = 2; % Guard sizes

    % Physics Constants
    dt = 1/fs; 
    df = fs/(M*N);
    Lfft_Sim = 2048;   

    % --- 1. TRANSMITTER SETUP (Local Logic) ---
    % Define masks LOCALLY to ensure 100% match
    [dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
    maskGuard = (abs(dd_L-lp)<=lTau) | (abs(dd_K-kp)<=2*kNu); 
    maskPilot = (dd_L==lp) & (dd_K==kp);
    maskData = ~(maskGuard | maskPilot);
    numData = nnz(maskData);

    % Local Helper to Build Grid
    function X = make_grid(data_syms)
        X = zeros(N, M);
        X(kp+1, lp+1) = Xp; % Pilot (MATLAB 1-based index)
        X(maskData) = data_syms;
    end

    % Verify Data Extraction (Safety Check)
    test_data = randn(numData, 1);
    X_test = make_grid(test_data);
    extracted = X_test(maskData);
    if norm(extracted - test_data) > 1e-10
        error('FATAL: Mask Mismatch! Cannot recover data from grid.');
    else
        fprintf('PASS: Grid Construction & Extraction aligned.\n');
    end

    % Generate Reference Grid
    rng(42);
    DataRef = qammod(randi([0 3], numData, 1), 4, 'UnitAveragePower', true);
    X_ref = make_grid(DataRef);

    % =========================================================================
    % FIGURE 5.1: Transmitted Delay-Doppler Grid
    % =========================================================================
    fprintf('Generating Fig 5.1 (Tx Grid)...\n');
    f1 = figure('Name', 'Fig 5.1: Transmitted Grid', 'Color', 'w', 'Position', [100 100 500 400]);
    imagesc(0:M-1, 0:N-1, abs(X_ref)); 
    axis xy; colormap(parula); colorbar;
    title('5.1 Transmitted Delay-Doppler Grid');
    xlabel('Delay Bin (l)'); ylabel('Doppler Bin (k)');
    clim([0 Xp]);

    % =========================================================================
    % FIGURE 5.3: Integer Channel Baseline
    % =========================================================================
    fprintf('Generating Fig 5.3 (Integer Baseline)...\n');
    tx_int = otfs_modulate(X_ref, N, M);
    if Mcp > 0, tx_int = [tx_int(end-Mcp+1:end); tx_int]; end

    params_int = struct('M',M, 'N',N, 'Mcp',Mcp, 'fs',fs);
    params_int.targets = struct('li',5, 'ki',4, 'gain',0.8); 

    rx_int = synth_channel_from_targets(tx_int, params_int);
    Y_int = otfs_demodulate(rx_int, Mcp, N, M);

    f3 = figure('Name', 'Fig 5.3: Integer Channel', 'Color', 'w', 'Position', [600 100 500 400]);
    imagesc(0:M-1, 0:N-1, abs(Y_int)); 
    axis xy; colormap(jet); colorbar;
    title('5.3 Demodulated Map (Integer Synthesized)');
    xlabel('Delay Bin'); ylabel('Doppler Bin');
    clim([0 Xp]);

    % =========================================================================
    % FIGURE 5.6: Fractional Alignment (Physics)
    % =========================================================================
    fprintf('Generating Fig 5.6 (Physics Alignment)...\n');
    t_frac = struct('gain', 0.8, 'tau', 5.3*dt, 'nu', 4.2*df);
    params_phy = params_int;
    params_phy.targets = t_frac;
    params_phy.Lfft = Lfft_Sim;

    tx_phy = otfs_modulate(X_ref, N, M);
    if Mcp > 0, tx_phy = [tx_phy(end-Mcp+1:end); tx_phy]; end 

    rx_phy = apply_channel(tx_phy, t_frac, params_phy);
    Y_phy = otfs_demodulate(rx_phy, Mcp, N, M);

    f6 = figure('Name', 'Fig 5.6: Physics Alignment', 'Color', 'w', 'Position', [100 550 500 400]);
    imagesc(0:M-1, 0:N-1, abs(Y_phy)); 
    axis xy; colormap(jet); colorbar;
    title('5.6 Demodulated Map (Corrected Alignment)');
    xlabel('Delay Bin'); ylabel('Doppler Bin');
    clim([0 Xp]);

    % =========================================================================
    % FIGURE 5.2: High-Resolution RVM (Sensing)
    % =========================================================================
    fprintf('Generating Fig 5.2 (High-Res RVM)...\n');
    Y_dual = fftshift(fft2(Y_phy));
    pad = 4;
    [Ny, Mx] = size(Y_phy);
    Y_padded = zeros(pad*Ny, pad*Mx);
    r_s = floor((pad*Ny - Ny)/2) + 1;
    c_s = floor((pad*Mx - Mx)/2) + 1;
    Y_padded(r_s:r_s+Ny-1, c_s:c_s+Mx-1) = Y_dual;
    RVM = abs(ifft2(ifftshift(Y_padded)));
    RVM = RVM / max(RVM(:));

    f2 = figure('Name', 'Fig 5.2: High-Res RVM', 'Color', 'w', 'Position', [600 550 500 400]);
    d_ax = linspace(0, M-1, size(RVM,2));
    k_ax = linspace(0, N-1, size(RVM,1));
    imagesc(d_ax, k_ax, RVM); 
    axis xy; colormap(jet); colorbar;
    title('5.2 High-Resolution RVM (Interpolated)');
    xlabel('Delay Bin'); ylabel('Doppler Bin');
    hold on; plot(lp + 5.3, kp + 4.2, 'wx', 'MarkerSize',12, 'LineWidth',2);

    % =========================================================================
    % FIGURE 5.4: Equalized Grid (MMSE)
    % =========================================================================
    fprintf('Generating Fig 5.4 (Equalization)...\n');
    % Build Heff
    Heff = build_Heff(N, M, params_phy, []); 

    noise_var_vis = 1e-4;
    Xhat_vec = mmse_equalize(Y_phy, Heff, noise_var_vis);
    Xhat = reshape(Xhat_vec, N, M);

    % Normalize
    pe = Xhat(kp+1, lp+1);
    if abs(pe)>0.1, Xhat = Xhat * (Xp / pe); end

    f4 = figure('Name', 'Fig 5.4: Equalized Grid', 'Color', 'w', 'Position', [1100 100 500 400]);
    imagesc(0:M-1, 0:N-1, abs(Xhat)); 
    axis xy; colormap(parula); colorbar;
    title('5.4 Equalized Grid (MMSE)');
    xlabel('Delay Bin'); ylabel('Doppler Bin');
    clim([0 Xp]);

    % =========================================================================
    % FIGURE 5.5: BER vs SNR (Rigorous Monte Carlo)
    % =========================================================================
    fprintf('Generating Fig 5.5 (Full Monte Carlo BER)...\n');
    SNR_dB = 0:4:20;
    BER = zeros(size(SNR_dB));
    MaxFrames = 500;  
    MinErrors = 200;

    sig_pwr = mean(abs(rx_phy).^2);

    for i = 1:length(SNR_dB)
        snr = SNR_dB(i);
        errs = 0; bits = 0; frame_cnt = 0;
        n_pwr = sig_pwr / (10^(snr/10));
        
        while (errs < MinErrors && frame_cnt < MaxFrames)
            frame_cnt = frame_cnt + 1;
            
            % 1. FRESH Data (Using Local Builder)
            DataNew = qammod(randi([0 3], numData, 1), 4, 'UnitAveragePower', true);
            X_MC = make_grid(DataNew);
            
            % 2. FRESH Transmit
            tx_MC = otfs_modulate(X_MC, N, M);
            if Mcp > 0, tx_MC = [tx_MC(end-Mcp+1:end); tx_MC]; end 
            
            % 3. Channel & Noise
            rx_MC = apply_channel(tx_MC, t_frac, params_phy);
            noise = sqrt(n_pwr/2) * (randn(size(rx_MC)) + 1j*randn(size(rx_MC)));
            y_MC = otfs_demodulate(rx_MC + noise, Mcp, N, M);
            
            % 4. Equalize
            x_eq_vec = mmse_equalize(y_MC, Heff, n_pwr);
            x_eq = reshape(x_eq_vec, N, M);
            
            % Normalize
            pe = x_eq(kp+1, lp+1);
            if abs(pe) > 0.01, x_eq = x_eq * (Xp/pe); end
            
            % 5. Error Count (Using Same Mask!)
            rx_syms = x_eq(maskData);
            tx_syms = X_MC(maskData);
            [~, ratio] = biterr(qamdemod(tx_syms,4), qamdemod(rx_syms,4), 2);
            
            errs = errs + (ratio * 2 * numel(tx_syms));
            bits = bits + (2 * numel(tx_syms));
            
            if snr >= 16 && errs == 0 && frame_cnt >= 20, break; end 
        end
        BER(i) = errs / bits;
        fprintf('SNR %2d dB: Frames=%3d, BER=%.4e\n', snr, frame_cnt, BER(i));
    end

    f5 = figure('Name', 'Fig 5.5: BER vs SNR', 'Color', 'w', 'Position', [1100 550 500 400]);
    semilogy(SNR_dB, BER, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    grid on; title('5.5 BER vs SNR');
    xlabel('SNR (dB)'); ylabel('Bit Error Rate');
    ylim([1e-6 1]);

    % =========================================================================
    % FIGURE 5.7: Constellation Diagram
    % =========================================================================
    fprintf('Generating Fig 5.7 (Constellation)...\n');
    maskConst = (abs(X_MC) > 0.1) & (abs(X_MC) < 9);
    tx_s = X_MC(maskConst);
    rx_s = x_eq(maskConst); 

    f7 = figure('Name', 'Fig 5.7: Constellation', 'Color', 'w', 'Position', [100 100 500 500]);
    plot(real(rx_s), imag(rx_s), 'r.', 'MarkerSize', 6); hold on;
    plot(real(tx_s), imag(tx_s), 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);
    grid on; axis equal; title(sprintf('5.7 Constellation (SNR %d dB)', snr));
    legend('Received', 'Transmitted');
    
    fprintf('Done. All figures generated.\n');
end