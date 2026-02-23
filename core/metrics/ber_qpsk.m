function [ber, total_errors, total_bits] = ber_qpsk(tx_syms, rx_syms)
% BER_QPSK Shared BER calculator.
rx_bits = qamdemod(rx_syms, 4, 'UnitAveragePower', true);
tx_bits = qamdemod(tx_syms, 4, 'UnitAveragePower', true);
[~, ber_ratio] = biterr(tx_bits, rx_bits, 2);
total_bits = 2 * numel(tx_bits);
total_errors = ber_ratio * total_bits;
ber = total_errors / max(total_bits,1);
end
