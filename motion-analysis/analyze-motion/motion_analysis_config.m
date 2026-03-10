function cfg = motion_analysis_config()
% MOTION_ANALYSIS_CONFIG  Default configuration for analyze_motion_and_QC.
%
%   cfg = motion_analysis_config() returns a struct containing all tunable
%   parameters for the motion analysis pipeline. Modify fields as needed
%   before passing the struct to analyze_motion_and_QC:
%
%   EXAMPLE:
%     cfg = motion_analysis_config();
%     cfg.pixel_size_um = 0.5;        % your pixel size in µm
%     cfg.fps = 2;                     % your frame rate in Hz
%     cfg.file_pattern = 'fly*.avi';   % match your file naming
%     analyze_motion_and_QC('/path/to/data', cfg);
%
%   See also: analyze_motion_and_QC

% =========================================================================
%  IMAGING PARAMETERS
% =========================================================================
% Pixel size in micrometers. This is used for all conversions from pixel
% shifts to physical units. Determine this from your optical system
% (e.g., sensor pixel pitch / magnification).
cfg.pixel_size_um   = 266 / 1024;   % µm/pixel (default: 50x with 200 mm tube lens, 2x2 binning)

% Acquisition frame rate in Hz. Must match your recording settings.
cfg.fps             = 1;            % Hz (default: 1 Hz for GFP movement recordings)

% Glob pattern for video file discovery. Only files matching this pattern
% in each data folder will be processed.
cfg.file_pattern    = '*.avi';      % e.g., '*.avi', 'S*.avi', 'fly*.tif'

% =========================================================================
%  NORMCORRE REGISTRATION PARAMETERS
% =========================================================================
% How the maximum search window is determined:
%   'auto'             - 30% of FOV, with rescue expansion up to 50%
%   'auto_conservative'- Fixed 50 µm window
%   'auto_loose'       - 50% of FOV from start
%   'fixed_um'         - User-specified in µm (requires max_shift_value)
%   'fixed_px'         - User-specified in pixels (requires max_shift_value)
cfg.max_shift_mode  = 'auto';

% Required for 'fixed_um' or 'fixed_px' modes. Scalar or [y, x] vector.
cfg.max_shift_value = [];

% Upsampling factor for subpixel registration accuracy. Higher values
% give finer subpixel resolution but increase computation time.
cfg.us_fac          = 50;

% Template update bin width in frames. The running template is updated
% by averaging this many consecutive frames.
cfg.bin_width       = 600;

% Number of frames used to build the initial reference template.
cfg.init_batch      = 600;

% Number of frames processed per chunk during streaming registration.
cfg.chunk_frames    = 600;

% RAM budget per processing chunk in GB. Chunk size is automatically
% reduced if the estimated memory exceeds this limit.
cfg.target_gb       = 2.0;

% =========================================================================
%  QUALITY CONTROL — HARD THRESHOLDS (gating)
% =========================================================================
% These thresholds determine the QC_PASS flag. Recordings that fail any
% hard threshold are flagged as unreliable for downstream analysis.

% Maximum fraction of frames where the estimated shift sits at the
% search window boundary ("pegged"). Values above this indicate that
% the search window was too small to track the actual motion.
cfg.qc_max_pegged_frac   = 0.01;   % fraction (1%)

% Maximum duration (seconds) of any single contiguous run of pegged
% frames. Even if the overall pegged fraction is low, a long
% contiguous run indicates a period of unreliable tracking.
cfg.qc_peg_run_min_s     = 10;     % seconds

% Maximum fraction of frames with NaN or non-finite shift values.
cfg.qc_nan_frac_max      = 0.01;   % fraction (1%)

% =========================================================================
%  QUALITY CONTROL — SOFT THRESHOLDS (descriptive, non-gating)
% =========================================================================
% These thresholds are reported in the output tables for informational
% purposes but do NOT affect the QC_PASS flag.

% Benchmark for "low motion" — median step size below this value
% is considered well-immobilized (soft criterion only).
cfg.qc_good_post_med_um  = 0.20;   % µm

% Spike detection: frames with step rate above this threshold, sustained
% for at least 2 seconds, are counted as spikes.
cfg.qc_spike_thr_um_s    = 0.50;   % µm/s

% Benchmark for spike fraction (informational).
cfg.qc_max_spike_frac    = 0.01;   % fraction

% Implausible jump detection threshold and fraction (informational).
cfg.qc_imp_jump_um_s     = 10.0;   % µm/s
cfg.qc_imp_jump_frac_max = 0.01;   % fraction

% =========================================================================
%  PROCESSING BEHAVIOR
% =========================================================================
% If true, skip files that already have a corresponding _shifts.csv.
% Set to false to force re-processing.
cfg.skip_if_present    = true;

% If true, include files containing 'cropped' in the filename.
cfg.include_cropped    = false;

% If true, use the Unicode character µ in plot axis labels.
% Set to false if your system does not render Unicode correctly.
cfg.use_unicode_labels = true;

% Random number generator seed for reproducibility of any stochastic
% operations (currently used for RNG state tracking in provenance).
cfg.rng_seed           = 0;

end