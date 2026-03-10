function cfg = neuro_motion_config()
% NEURO_MOTION_CONFIG  Default configuration for analyze_motion_and_QC_neuro.
%
%   cfg = neuro_motion_config() returns a struct containing all tunable
%   parameters for the neurological motion analysis pipeline. Modify fields
%   as needed before passing the struct to analyze_motion_and_QC_neuro:
%
%   EXAMPLE:
%     cfg = neuro_motion_config();
%     cfg.pixel_size_um = 0.5;        % your pixel size in µm
%     cfg.fps = 30;                    % your frame rate in Hz
%     cfg.file_pattern = 'fly*.avi';   % match your file naming
%     analyze_motion_and_QC_neuro('/path/to/data', cfg);
%
%   See also: analyze_motion_and_QC_neuro

% =========================================================================
%  IMAGING PARAMETERS
% =========================================================================
% Pixel size in micrometers. This is used for all conversions from pixel
% shifts to physical units. Determine this from your optical system
% (e.g., sensor pixel pitch / magnification).
cfg.pixel_size_um   = 266 / 1024;   % µm/pixel (default: 50x with 200 mm tube lens, 2x2 binning)

% Nominal acquisition frame rate in Hz. Must match your recording settings.
% At 33.33 Hz, 10,000 frames = 300 seconds = 5.00 minutes exactly.
cfg.fps             = 1000 / 30;    % Hz (default: 33.33 Hz for calcium recordings)

% Override for AVI-reported frame rate. AVI containers typically truncate
% fps to integer (33 Hz instead of 33.33 Hz). Set to [] to trust
% VideoReader's reported frame rate instead.
cfg.fps_override    = 1000 / 30;    % Hz (default: 33.33 Hz; set [] to use VideoReader)

% Glob pattern for video file discovery. Only files matching this pattern
% in each data folder will be processed.
cfg.file_pattern    = 'S*.avi';     % e.g., '*.avi', 'S*.avi', 'fly*.tif'

% =========================================================================
%  NORMCORRE REGISTRATION PARAMETERS
% =========================================================================
% How the maximum search window is determined:
%   'auto'             - 40% of min(FOV) (handles saccades; default)
%   'auto_conservative'- 10% of min(FOV)
%   'auto_loose'       - 45% of min(FOV)
%   'auto_extreme'     - 50% of min(FOV) (last resort)
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
cfg.bin_width       = 240;

% Number of frames used to build the initial reference template.
cfg.init_batch      = 300;

% Number of frames processed per chunk during streaming registration.
cfg.chunk_frames    = 400;

% RAM budget per processing chunk in GB. Chunk size is automatically
% reduced if the estimated memory exceeds this limit.
cfg.target_gb       = 2.0;

% =========================================================================
%  QUALITY CONTROL — MOVEMENT_VALID THRESHOLDS (tracking integrity)
% =========================================================================
% These thresholds determine the MOVEMENT_VALID flag. Recordings that
% fail any hard threshold have unreliable motion tracking and cannot be
% used for immobilization comparison.

% Maximum fraction of frames at the search window boundary ("pegged").
% Values above this indicate the search window was too small.
cfg.qc_max_pegged_frac   = 0.20;   % fraction (20%)

% Maximum fraction of frames with NaN or non-finite shift values.
cfg.qc_nan_frac_max      = 0.05;   % fraction (5%)

% Flatline and unstable_tracking always fail MOVEMENT_VALID (no threshold).

% =========================================================================
%  QUALITY CONTROL — ROI_ELIGIBLE THRESHOLDS (neural pipeline gating)
% =========================================================================
% ROI_ELIGIBLE requires MOVEMENT_VALID = true PLUS boundary checks.
% If the brain exits the FOV, motion tracking is still valid for
% immobilization comparison but neural data would be incomplete.

% FOV excursion fraction for warning (informational only).
cfg.boundary_warning_frac = 0.25;   % fraction (25% of FOV)

% FOV excursion fraction for ROI ineligibility. If total excursion in
% either axis exceeds this fraction of the FOV, ROI_ELIGIBLE = false.
cfg.boundary_fail_frac    = 0.40;   % fraction (40% of FOV)

% =========================================================================
%  QUALITY CONTROL — SOFT THRESHOLDS (descriptive, non-gating)
% =========================================================================
% These are reported for informational purposes but do NOT affect either
% MOVEMENT_VALID or ROI_ELIGIBLE flags.

% Benchmark for "low motion" — median step size below this value
% is considered well-immobilized (informational only).
cfg.qc_good_post_med_um  = 0.30;   % µm

% Spike detection: step rate above this threshold, sustained for a
% minimum number of consecutive frames, is counted as a spike.
cfg.qc_spike_thr_um_s    = 3.0;    % µm/s

% Minimum consecutive frames for spike detection. If empty, defaults
% to ~180 ms worth of frames (max(6, round(0.18 * fps))).
cfg.qc_spike_min_consec_fr = [];

% Benchmark for spike fraction (informational).
cfg.qc_max_spike_frac    = 0.01;   % fraction

% =========================================================================
%  REPOSITIONING DETECTION
% =========================================================================
% Detects sustained large-step events with a baseline shift, indicating
% the sample was physically repositioned during recording.

% Step size threshold for candidate repositioning frames.
cfg.reposition_thr_um          = 2.0;   % µm

% Minimum duration of sustained above-threshold motion to qualify.
cfg.reposition_min_sustain_s   = 2.0;   % seconds

% Minimum baseline position shift (pre vs post) to confirm repositioning.
cfg.reposition_baseline_shift_um = 2.0; % µm

% =========================================================================
%  SETTLE POLICY
% =========================================================================
% Post-settle metrics exclude early frames to allow the preparation to
% stabilize. The effective settle period is:
%   min(base_settle_min, max(min_settle_min, frac_settle * duration_min))

cfg.settle_base_min   = 10;    % minutes (maximum settle period)
cfg.settle_min_min    = 0.5;   % minutes (minimum settle period)
cfg.settle_frac       = 0.20;  % fraction of recording duration

% =========================================================================
%  PROCESSING BEHAVIOR
% =========================================================================
% If true, skip files that already have a corresponding _shifts.csv.
% Set to false to force re-processing.
cfg.skip_if_present    = true;

% If true, use the Unicode character µ in plot axis labels.
% Set to false if your system does not render Unicode correctly.
cfg.use_unicode_labels = true;

% Random number generator seed for reproducibility of any stochastic
% operations (currently used for RNG state tracking in provenance).
cfg.rng_seed           = 0;

end