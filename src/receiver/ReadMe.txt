include / ensure the following fields inside waveform_package.mat (exact names used by the code):

params.M and params.N — OTFS grid sizes (or ensure X is present so Person B can infer).

params.Mcp — cyclic prefix length (0 if none).

params.Lfft — FFT length used for channel synthesis (helpful).

params.fs (optional) — sampling rate (not required by current code but useful).

params.targets — array of target structs. Each target should have:

li (integer delay bin), ki (integer Doppler bin), and gain (complex).

If different names used, include synonyms (delay_idx, doppler_idx).

(Strongly recommended) params.txch_rx_map — the precomputed effective linear map (size MN × MN). If you provide this matrix, build_Heff will load it directly and skip expensive simulation.

(Optional) apply_channel(tx, targets, params) function on MATLAB path — if Person A provides it, build_Heff will call it to build Heff by sending basis vectors through the true channel (accurate but slow).

X — the DD grid (pilot+guard+data).

tx_signal — the time-domain transmit signal (column vector).

params.noise_var (optional) — noise variance used by MMSE.

If params.txch_rx_map is present, Person B will use it immediately (fast). If not, the pipeline falls back to either an apply_channel function or a simple synthesized channel from params.targets.

TROUBLESHOOTING / TIPS

Heff construction via apply_channel or per-column simulation is slow (MN linear solves, may be heavy for 32×32 = 1024 columns). If you want speed, Person A should provide params.txch_rx_map.

If params.targets contains non-integer fractional delays/dopplers, the simple synthesizer will approximate fractional effects poorly — providing apply_channel or txch_rx_map is best.

If your OTFS modulation uses a different ISFFT/SFFT convention (row/column ordering or scaling), adjust otfs_demodulate and otfs_modulate_local accordingly to match Person A’s transmitter.