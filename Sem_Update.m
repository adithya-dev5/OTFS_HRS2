% Sem_Update_fixed.m
% Corrected single-file MATLAB script for Semester-1 OTFS plots and BER baseline.
% Fixes:
%  - Replaced interp1('sinc') call with FIR sinc fractional-delay filter (portable)
%  - Replaced implicit expansion with bsxfun for constellation distance (portable)
%  - Made de2bi usage robust by explicit orientations
%  - Minor robustness fixes for indexing and shapes
% Save as Sem_Update_fixed.m and run in MATLAB.

clear; close all; clc;

%% ----------------- Parameters (edit for your setup) --------------------
M = 32;                
N = 32;                
lp = floor(M/2)+1;     
kp = floor(N/2)+1;     
lTau = 2;              
kNu  = 2;              
QAMOrder = 4;          
numBits = 1e4;         
fc = 4e9;              
c = 3e8;               
v_ped = 1.0;           
lambda = c/fc;         
nu_Hz = 2*v_ped/lambda;

fractional_doppler = 0.05;  
fractional_delay   = 0.3;  

SNR_dB_for_heatmap = 30;     

CP_len = 8;                  

%% ----------------- Build DD grid with Pilot + Guard + Random Data -------------
Xdd = zeros(M, N);    

Xp = 5+0j; 
Xdd(lp, kp) = Xp;

rng(0); 
symset = pskmod(0:QAMOrder-1, QAMOrder, pi/QAMOrder); 

for r = 1:M
    for c = 1:N
        if (r==lp && c==kp), continue; end
        Xdd(r,c) = symset(randi(length(symset)));
    end
end

figure('Name','1. DD grid (pilot+guard+data)','NumberTitle','off');
imagesc(1:M,1:N,abs(Xdd).'); axis xy; axis tight;
xlabel('Delay bin l'); ylabel('Doppler bin k');
title('DD grid magnitude |X_{DD}(l,k)| (pilot at center)');
colorbar;
hold on;
plot(lp, kp, 'r+', 'MarkerSize',12, 'LineWidth',2);

%% ----------------- Simplified OTFS helpers (TF/DD transforms) ------------
ISFFT = @(Xdd) ifft(Xdd, [], 2);  
SFFT  = @(Xtf) fft(Xtf, [], 2);   

ofdm_mod = @(Xtf) ifft(Xtf, [], 1);    
ofdm_demod = @(td_blocks) fft(td_blocks, [], 1); 

make_tx_signal = @(td_blocks) reshape([td_blocks(end-CP_len+1:end, :); td_blocks], [], 1);

extract_td_blocks = @(rx_vec, Mloc, Nloc) reshape(rx_vec, Mloc+CP_len, Nloc);
strip_cp = @(blocks) blocks(CP_len+1:end, :);

%% ----------------- Transmit through simple fractional-delay & Doppler channel ----------
Xtf = ISFFT(Xdd);
td_blocks = ofdm_mod(Xtf);
tx_signal_vec = make_tx_signal(td_blocks);

l_int = 3;                   
frac = fractional_delay;     
k_int = 4;                   
frac_dop = fractional_doppler;

total_samples = length(tx_signal_vec);
t_idx = (0:total_samples-1).';

doppler_phase = exp(1j*2*pi*frac_dop * (t_idx/total_samples) * k_int * N);

int_delay_samples = l_int;
tx_delayed = [zeros(int_delay_samples,1); tx_signal_vec(1:end-int_delay_samples)];

L = 8; 
n = -L:L;
h = sinc(n - frac);
h = h ./ sum(h);
tx_frac_delay = conv(tx_delayed, h, 'same');

rx_signal_clean = tx_frac_delay .* doppler_phase;

sigma = 10^(-SNR_dB_for_heatmap/20);
rx_signal = rx_signal_clean + sigma*(randn(size(rx_signal_clean)) + 1j*randn(size(rx_signal_clean)))/sqrt(2);

%% ----------------- Receiver: OFDM demod -> SFFT -> DD domain -------------
rx_blocks = extract_td_blocks(rx_signal, M, N);
rx_nocp = strip_cp(rx_blocks);

RX_tf = ofdm_demod(rx_nocp);
Ydd = SFFT(RX_tf);
Ydd = Ydd / max(abs(Ydd(:)));

figure('Name','2. DD-domain received magnitude |Y(l,k)|','NumberTitle','off');
imagesc(1:M, 1:N, abs(Ydd).'); axis xy; axis tight;
xlabel('Delay bin l'); ylabel('Doppler bin k');
title('Received DD-domain magnitude |Y(l,k)|');
colorbar;

figure('Name','3. Range-Velocity heatmap |Y|^2','NumberTitle','off');
imagesc(1:M,1:N, abs(Ydd).'.^2 ); axis xy; axis tight;
xlabel('Delay bin l'); ylabel('Doppler bin k');
title('Range-Velocity heatmap R(l,k)=|Y|^2');
colorbar;

[~, idx_lin] = max(abs(Ydd(:)));
[delay_idx, dop_idx] = ind2sub(size(Ydd), idx_lin);
fprintf('Strongest bin at delay=%d, doppler=%d\n', delay_idx, dop_idx);

figure('Name','4. Doppler slice','NumberTitle','off');
plot(1:N, abs(Ydd(delay_idx, :)), '-o','LineWidth',1.5);
xlabel('Doppler bin k'); ylabel('|Y(l_{max}, k)|');
title(sprintf('Doppler slice at delay bin l=%d', delay_idx));
grid on;

%% ----------------- BER vs SNR baseline (AWGN) -------------------------
SNRdB_list = 0:2:20;
BER = zeros(size(SNRdB_list));
numFramesBER = 60;

sym_map = pskmod(0:QAMOrder-1, QAMOrder, pi/QAMOrder); 
numSym = M * N;

for si = 1:length(SNRdB_list)
    snrdB = SNRdB_list(si);
    totalErr = 0; totalBits = 0;
    for f = 1:numFramesBER

        txIdx = randi([0 QAMOrder-1], numSym,1);
        txSymbols = sym_map(txIdx+1);            
        txSymbols = reshape(txSymbols, M, N);

        Xtf_tx = ISFFT(txSymbols);
        td_blocks_tx = ofdm_mod(Xtf_tx);
        tx_vec = make_tx_signal(td_blocks_tx);

        noise_sigma = 10^(-snrdB/20);
        rx_vec_awgn = tx_vec + noise_sigma*(randn(size(tx_vec))+1j*randn(size(tx_vec)))/sqrt(2);

        rx_blocks_tmp = extract_td_blocks(rx_vec_awgn, M, N);
        rx_nocp_tmp = strip_cp(rx_blocks_tmp);
        RX_tf_tmp = ofdm_demod(rx_nocp_tmp);
        Ydd_rx = SFFT(RX_tf_tmp);

        rxSymbols_vec = reshape(Ydd_rx, [], 1);

        D = bsxfun(@minus, rxSymbols_vec, sym_map);
        dists = abs(D).^2;
        [~, idxHat] = min(dists, [], 2);
        idxHat = idxHat - 1;

        txBits = de2bi(txIdx(:), log2(QAMOrder), 'left-msb');
        rxBits = de2bi(idxHat(:), log2(QAMOrder), 'left-msb');
        bitErrs = sum(sum(txBits ~= rxBits));

        totalErr = totalErr + bitErrs;
        totalBits = totalBits + numSym * log2(QAMOrder);
    end
    BER(si) = totalErr / totalBits;
end

figure('Name','5. BER vs SNR (AWGN baseline)','NumberTitle','off');
semilogy(SNRdB_list, BER, '-s','LineWidth',1.6);
grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('OTFS baseline BER vs SNR (AWGN)');
