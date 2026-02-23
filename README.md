# OTFS ISAC Research Architecture (Hardened)

## Structure
- `config/`: canonical configuration (`default_config.m`)
- `core/tx|channel|rx|equalization|sensing|metrics`: reusable logic modules
- `models/otfs`, `models/ofdm`: model-level organization
- `experiments/`: non-plotting experiment runners
- `analysis/`: generated outputs and post-processing artifacts
- `tests/`: unit test scripts

## Run experiments
1. `addpath('config'); addpath('core/tx'); addpath('core/channel'); addpath('core/rx'); addpath('core/equalization'); addpath('core/sensing'); addpath('core/metrics'); addpath('core/utils'); addpath('experiments');`
2. `cfg = default_config();`
3. `r1 = run_ber_vs_snr(cfg);`
4. `r2 = run_rmse_vs_snr(cfg);`
5. `r3 = run_pd_vs_ber(cfg);`

If you forget to define `cfg`, you can call runners with defaults:
- `r1 = run_ber_vs_snr();`
- `r2 = run_rmse_vs_snr();`
- `r3 = run_pd_vs_ber();`

All runners save structured outputs to `cfg.experiments.output_dir` and return result structs.

## Tests
- `tests/test_dd_to_physical_mapping.m`
- `tests/test_channel_scaling.m`
- `tests/test_rvm_peak_location.m`
- `tests/run_all_tests.m` (from repo root: `run('tests/run_all_tests.m')` or `tests.run_all_tests` if package-scoped)

## Numerical invariance harness
- `core/metrics/validate_numerical_invariance.m` compares pre-refactor (`src/*`) and post-refactor (`core/*`) outputs for DD grid, channel output, RVM, BER, and RMSE under identical seed/config and halts if any deviation exceeds `1e-10`.
