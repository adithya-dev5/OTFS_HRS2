# Waveform & Channel utilities (NITC OTFS project)

Files (in this folder)
- waveform_and_channel.m      : top-level builder (returns X, tx, params, targets)
- build_dd_grid.m             : construct N x M DD grid with pilot/guard/data
- otfs_modulate.m             : ISFFT/IFTT chain (IFFT across rows, reshape), normalized tx
- apply_channel.m             : fractional-delay (freq-domain) + time-domain Doppler channel
- add_cp.m                    : prepend cyclic prefix (Mcp samples)
- save_tx_for_receiver.m      : helper to store MAT reference for Person B
- test_build_dd_grid.m        : unit test reproducing DD_Grid figure (masks & counts)
- test_tx_plot.m              : compute tx and plot abs(tx) (compare to Receiver.txt)
- test_channel_shift_test.m   : sanity-check applying single-target channel and verify DD peak
- package_and_save.m          : builds standard example and saves waveform_package.mat

MAT variables stored in waveform_package.mat
- X           : N x M delay-Doppler grid (DD domain)
- tx_signal   : normalized transmit waveform (no CP)
- params      : struct with fields (fs, M, N, Mcp, Lfft, fvec, tvec)
- Lfft        : FFT length used (same as params.Lfft)
- fvec        : frequency vector used for fractional-delay operations
- tvec        : time vector for Doppler multiplication
- targets     : example target array (struct with fields tau, nu, gain)

How to run tests
1. Add this folder to MATLAB path.
2. Run `test_build_dd_grid.m` — verify PASS messages and the DD-grid figure.
3. Run `test_tx_plot.m` — verify tx plot matches project Receiver plot (visual).
4. Run `test_channel_shift_test.m` — inspect PASS/WARN and the DD magnitude if mismatch.
5. Run `package_and_save.m` to generate waveform_package.mat for Person B.

Notes & conventions
- Indices lp,kp are treated as 0-based internally; the code tolerates typical MATLAB 1-based inputs heuristically.
- apply_channel requires `params.fs` (sampling frequency) to compute tvec and fvec.
- Fractional delays are implemented by multiplying the zero-padded FFT of tx by exp(-j*2*pi*f*tau) and inverse transforming (see apply_channel.m).
- These tests follow the DD-grid plotting / mask logic shown in Midterm_GC05 (DD_Grid snippet). :contentReference[oaicite:4]{index=4}
