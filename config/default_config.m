function cfg = default_config()
% DEFAULT_CONFIG Canonical configuration for OTFS-ISAC experiments.
% Single source of truth for radar, waveform, axis, channel, and experiment control.
%
% All scripts and core functions should consume this cfg struct (or a derivative)
% instead of hardcoding parameters.

% --------------------------
% Reproducibility
% --------------------------
cfg.seed = 42;

% --------------------------
% Radar / waveform parameters
% --------------------------
cfg.fc = 77e9;           % Carrier frequency (Hz)
cfg.B = 20e6;            % Bandwidth (Hz)
cfg.M = 32;              % Delay bins (columns)
cfg.N = 32;              % Doppler bins (rows)
cfg.Mcp = 0;             % Cyclic prefix length (samples)
cfg.fs = 15.36e6;        % Sampling rate (Hz)

% Derived timing/frequency quantities (kept explicit for consistency)
cfg.Delta_f = cfg.B / cfg.M;           % Subcarrier spacing (Hz)
cfg.T = 1 / cfg.Delta_f;               % Symbol duration (s)
cfg.Ts = 1 / (cfg.M * cfg.Delta_f);    % Sample interval (s)

% --------------------------
% DD grid and guard geometry
% --------------------------
cfg.lp = floor(cfg.M/2) - 1;  % Pilot delay index (0-based)
cfg.kp = floor(cfg.N/2) - 1;  % Pilot doppler index (0-based)
cfg.lTau = 2;                 % Guard half width in delay
cfg.kNu = 2;                  % Guard half width in doppler (used as 2*kNu)
cfg.Xp = 10;                  % Pilot amplitude

% --------------------------
% Axis and convention control
% --------------------------
cfg.dd_axis = "doppler_row_delay_col"; % Mandatory convention

% --------------------------
% Channel selection and schema
% --------------------------
cfg.channel.mode = "synthetic_fractional"; % synthetic_integer|synthetic_fractional|backscatter
cfg.channel.Lfft = 2048;
cfg.channel.noise_snr_db = 30;

% Target schema (single source): struct with fields tau (s), nu (Hz), gain (complex)
cfg.targets = struct('tau', 5/cfg.fs, 'nu', 4*(cfg.fs/(cfg.M*cfg.N)), 'gain', 0.8);

% --------------------------
% Normalization flags
% --------------------------
cfg.normalize.tx = false;          % keep false to preserve current Heff linearity behavior
cfg.normalize.rx_pilot_rescale = true;

% --------------------------
% Equalization configuration
% --------------------------
cfg.equalization.method = "mmse";
cfg.equalization.solver = "direct"; % direct|iterative (iterative kept as optional stub)

% --------------------------
% Experiment defaults
% --------------------------
cfg.experiments.snr_db_range = 0:2:20;
cfg.experiments.max_frames = 500;
cfg.experiments.min_errors = 300;
cfg.experiments.output_dir = 'analysis/results';

% --------------------------
% Runtime / integrity controls
% --------------------------
cfg.enable_profiling = false;
cfg.git_commit_hash = '';
cfg.seed_initialized = false;

end
