% compare_tx_rx.m
% Compare transmitted and received signals after running run_integration.
% Assumes waveform_package.mat (from Person A) and receiver_results.mat exist.

clear; close all; clc;

% --- load files
wp = load('waveform_package.mat');
R  = load('receiver_results.mat');

% tx: try common names
if isfield(wp,'tx_signal'), tx = wp.tx_signal(:);
elseif isfield(wp,'tx'), tx = wp.tx(:);
else error('tx_signal not found in waveform_package.mat'); end

% rx: try from receiver results or build from rx saved var
if isfield(R,'rx'), rx = R.rx(:);
elseif isfield(R,'y'), rx = R.y(:);
elseif isfield(R,'Y') && isfield(wp,'params') && isfield(wp.params,'M') && isfield(wp.params,'N')
    % we may reconstruct time-domain rx used in run_integration by re-modulating Y_tf if saved
    % but typical case: receiver_results contains no time-domain rx; try to load file 'rx_signal' if exists
    if isfield(R,'rx_signal'), rx = R.rx_signal(:);
    else
        % fallback: try reconstruct rx from Xhat via OTFS modulation (approx)
        fprintf('No rx vector in receiver_results.mat; will try reconstructing approximate rx from Xhat if present.\n');
        if isfield(R,'Xhat')
        % --- FIX: Use actual otfs_modulate to reconstruct time domain ---
        Xhat = R.Xhat;
        % Ensure we have N and M
        if isfield(wp.params,'N'), N=wp.params.N; M=wp.params.M; else [N,M]=size(Xhat); end
        rx = otfs_modulate(Xhat, N, M);
        fprintf('Reconstructed rx from Xhat using otfs_modulate.\n');
        
        else
            error('No rx available in receiver_results.mat and no Xhat to reconstruct. Cannot proceed.');
        end
    end
end

% If CP was used, remove CP for direct comparison to tx (or add if needed)
params = wp.params;
if isfield(params,'Mcp'), Mcp = params.Mcp; else Mcp = 0; end

% Make sure tx and rx are similar lengths by truncating/padding
Ltx = numel(tx);
Lrx = numel(rx);
L = max(Ltx, Lrx);

txv = tx; rxv = rx;
if Ltx < L, txv = [txv; zeros(L-Ltx,1)]; end
if Lrx < L, rxv = [rxv; zeros(L-Lrx,1)]; end

% Time axis (samples)
n = (0:L-1).';

% --- 1) Time-domain plots (real and imag)
figure('Name','Time domain (real)'); 
subplot(2,1,1);
plot(n, real(txv), '-'); hold on; plot(n, real(rxv), '--'); hold off;
legend('tx (real)','rx (real)'); xlabel('sample'); title('Real parts: tx vs rx'); grid on;

subplot(2,1,2);
plot(n, imag(txv), '-'); hold on; plot(n, imag(rxv), '--'); hold off;
legend('tx (imag)','rx (imag)'); xlabel('sample'); title('Imag parts: tx vs rx'); grid on;

% --- 2) Envelope (magnitude) and instantaneous phase
figure('Name','Envelope & Phase');
subplot(2,1,1); plot(n, abs(txv),'b-'); hold on; plot(n, abs(rxv),'r.--'); hold off;
legend('|tx|','|rx|'); title('Envelope (magnitude)'); xlabel('sample'); grid on;
subplot(2,1,2); plot(n, unwrap(angle(txv)),'b-'); hold on; plot(n, unwrap(angle(rxv)),'r.--'); hold off;
legend('phase tx','phase rx'); title('Instantaneous phase'); xlabel('sample'); grid on;

% --- 3) Complex correlation (cross-correlation) -> delay estimate
[cc, lags] = xcorr(rxv, txv, 'coeff'); % normalized cross-correlation
[~, idxmax] = max(abs(cc));
delay_est = lags(idxmax);
fprintf('Estimated sample delay (rx relative to tx): %d samples (positive means rx delayed).\n', delay_est);

figure('Name','Cross-correlation');
plot(lags, abs(cc)); grid on; title('Cross-correlation magnitude (rx vs tx)');
xlabel('lag (samples)'); ylabel('abs(xcorr)'); hold on;
plot(delay_est, abs(cc(idxmax)), 'ro','MarkerSize',10); hold off;

% --- 4) MSE, SNR estimate and complex inner-product (gain/phase)
% Align rx to tx using delay estimate
if delay_est >= 0
    rx_al = rxv( (1+delay_est) : min(L, L + delay_est) );
    tx_al = txv(1 : numel(rx_al));
else
    % negative delay => rx leads tx
    dd = abs(delay_est);
    tx_al = txv( (1+dd) : min(L, L + dd) );
    rx_al = rxv(1 : numel(tx_al));
end

% ensure same length
Nal = min(numel(tx_al), numel(rx_al));
tx_al = tx_al(1:Nal); rx_al = rx_al(1:Nal);

% compute complex scale (least-squares) to align amplitude/phase: find alpha minimizing ||rx - alpha*tx||
alpha = (tx_al' * tx_al) \ (tx_al' * rx_al);   % alpha = (tx^H tx)^{-1} tx^H rx
rx_pred = alpha * tx_al;
resid = rx_al - rx_pred;
mse = mean(abs(resid).^2);
power_signal = mean(abs(tx_al).^2) * abs(alpha)^2;
est_snr = 10*log10(power_signal / mse);
fprintf('Complex scale alpha (amplitude+phase) = %.4f + %.4fj (abs=%.3f, angle=%.3f rad)\n', real(alpha), imag(alpha), abs(alpha), angle(alpha));
fprintf('MSE (after alignment & scaling) = %.4e. Estimated SNR ≈ %.2f dB\n', mse, est_snr);

% plot aligned signals over zoomed window
zoomRange = 1:min(200, Nal);
figure('Name','Aligned signals (zoom)');
plot(zoomRange, real(tx_al(zoomRange)),'b-'); hold on;
plot(zoomRange, real(rx_al(zoomRange)),'r--'); plot(zoomRange, real(rx_pred(zoomRange)),'k:'); hold off;
legend('tx','rx (aligned)','alpha*tx'); title('Aligned real parts (zoomed)'); grid on;

% --- 5) Frequency-domain (FFT) comparison / PSD
Nfft = 2^nextpow2(max(4096, Nal));
TXF = fftshift(fft(txv, Nfft));
RXF = fftshift(fft(rxv, Nfft));
freq = linspace(-0.5,0.5,Nfft); % normalized freq (cycles/sample)

figure('Name','Frequency-domain magnitude (normalized freq)');
plot(freq, 20*log10(abs(TXF)/max(abs(TXF)))); hold on;
plot(freq, 20*log10(abs(RXF)/max(abs(RXF)))); hold off;
legend('tx (dB normalized)','rx (dB normalized)'); xlabel('Normalized freq (cycles/sample)'); ylabel('dB'); grid on;
title('Normalized spectra');

% --- 6) Spectrograms to visualize Doppler / time-varying freq
figure('Name','Spectrograms (tx, rx)');
subplot(2,1,1);
spectrogram(txv, 256, 200, 256, 1, 'yaxis'); title('Spectrogram: tx');
subplot(2,1,2);
spectrogram(rxv, 256, 200, 256, 1, 'yaxis'); title('Spectrogram: rx');

% --- 7) Constellation comparison (if QPSK-like transmitted X exists)
if isfield(wp,'X') && isfield(R,'Xhat')
    Xtx = wp.X;
    Xrx = R.Xhat;
    % pick data positions (where Xtx nonzero)
    data_mask = (abs(Xtx) > 0);
    tx_sym = Xtx(data_mask);
    rx_sym = Xrx(data_mask);
    figure('Name','Constellation: tx vs rx'); 
    plot(real(tx_sym(:)), imag(tx_sym(:)), 'bo'); hold on;
    plot(real(rx_sym(:)), imag(rx_sym(:)), 'rx'); hold off;
    legend('tx','rx'); xlabel('Real'); ylabel('Imag'); title('DD-domain symbols: tx vs rx'); grid on;
end

% --- 8) Save summary numbers
summary.delay_samples = delay_est;
summary.alpha = alpha;
summary.mse = mse;
summary.estimated_snr_db = est_snr;
save('tx_rx_comparison_summary.mat', 'summary');

fprintf('Comparison complete. Summary saved to tx_rx_comparison_summary.mat\n');
