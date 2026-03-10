function cfg = compare_motion_config()
% COMPARE_MOTION_CONFIG  Default configuration for compare_motion_groups.
%
%   cfg = compare_motion_config() returns a struct containing all tunable
%   parameters. Modify fields as needed before passing to the comparison:
%
%   EXAMPLE:
%     cfg = compare_motion_config();
%     cfg.control_label = 'Saline';
%     cfg.experimental_label = 'Drug X';
%     cfg.time_window_min = [0 30];     % analyze first 30 min only
%     compare_motion_groups('/ctrl', '/exp', '/output', cfg);
%
%   See also: compare_motion_groups

% =========================================================================
%  GROUP LABELS (used in plots, tables, and filenames)
% =========================================================================
cfg.control_label      = 'Hydrogel';
cfg.experimental_label = 'Diethyl Ether + Hydrogel';

% =========================================================================
%  TIME WINDOW
% =========================================================================
% Restrict analysis to a specific time range (in minutes).
%   [0 Inf]  = full recording (default)
%   [0 30]   = first 30 minutes only
%   [10 60]  = minutes 10 through 60
% The window is applied to the per-frame shift CSVs using the time_min
% column. Frames outside the window are excluded before metric computation.
cfg.time_window_min    = [0 Inf];   % [start_min, end_min]

% =========================================================================
%  STATISTICAL PARAMETERS
% =========================================================================
% Number of bootstrap resamples for confidence intervals and permutation
% tests. 10,000 is standard for publication; reduce for faster iteration.
cfg.n_bootstrap        = 10000;

% Random number generator seed for reproducibility.
cfg.rng_seed           = 42;

% Minimum sample size before warnings are issued.
cfg.min_sample_size    = 5;

% =========================================================================
%  EXPORT SETTINGS
% =========================================================================
% Resolution for PNG figure export.
cfg.export_dpi         = 300;

% If true, also export PDF (vector) versions of all figures.
cfg.export_vector      = true;

% If true, use robust Y-axis limits (clip outliers beyond 1.5*IQR).
cfg.robust_ylim        = true;

end