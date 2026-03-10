function cfg = roi_refinement_config()
% ROI_REFINEMENT_CONFIG  Default configuration for ROI_REFINEMENT_v3_2.
%
%   cfg = roi_refinement_config() returns a struct containing all tunable
%   parameters for the ROI refinement pipeline. Modify fields as needed
%   before passing to ROI_REFINEMENT_v3_2:
%
%   EXAMPLE:
%     cfg = roi_refinement_config();
%     cfg.conditions(1).name = 'Control';
%     cfg.conditions(1).roi_analysis_dir = '/path/to/roi_analysis';
%     cfg.conditions(1).output_dir = '/path/to/refined';
%     cfg.summary_output_dir = '/path/to/summary';
%     ROI_REFINEMENT_v3_2(cfg);
%
%   If conditions are not specified, a folder picker opens.
%
%   See also: ROI_REFINEMENT_v3_2

% =========================================================================
%  DATA LOCATIONS
% =========================================================================
% Each condition needs:
%   .name             - Label (e.g., 'Control', 'Experimental')
%   .roi_analysis_dir - Folder containing *_rois.mat files from pipeline
%   .output_dir       - Folder for refined output
%
% If left empty, the pipeline will prompt for folder selection.
cfg.conditions = struct('name', {}, 'roi_analysis_dir', {}, 'output_dir', {});

% Master output for cross-condition summaries. If empty, uses the
% parent directory of the first condition's output_dir.
cfg.summary_output_dir = '';

% =========================================================================
%  IMAGING PARAMETERS (must match pipeline)
% =========================================================================
cfg.frame_rate    = 33.33;   % Hz
cfg.pixel_size_um = 0.260;   % µm/pixel

% =========================================================================
%  SHAPE CRITERIA (v3.0)
% =========================================================================
% Edge contact: Remove ROIs touching brain mask edge
cfg.max_edge_contact_fraction = 0.20;  % Remove if >20% of perimeter on edge

% Compactness: 4*pi*area/perimeter^2 (circle=1, complex shapes <1)
cfg.min_compactness = 0.20;  % Remove very irregular shapes

% Bounding box extent: area / bounding_box_area
cfg.min_extent = 0.35;  % Remove stringy ROIs
cfg.max_extent = 0.95;  % Remove perfectly rectangular ROIs (artifacts)

% Aspect ratio of bounding box
cfg.max_aspect_ratio = 3.0;

% Solidity
cfg.min_solidity = 0.55;

% Eccentricity
cfg.max_eccentricity = 0.93;

% =========================================================================
%  TRACE / KINETICS CRITERIA (v3.0)
% =========================================================================
% Negative transients
cfg.max_negative_zscore = -3.0;  % Remove if dips > 3 std below baseline

% Rise/decay ratio (GCaMP8f: fast rise ~50ms, slow decay ~200-500ms)
cfg.max_rise_decay_ratio = 3.0;

% Step artifacts: instantaneous jumps
cfg.max_single_frame_jump = 0.5;  % dF/F change in one frame

% Sustained plateaus
cfg.max_plateau_duration_s = 20.0;
cfg.plateau_threshold = 0.7;  % Fraction of peak that counts as "elevated"

% Decay time constant range (GCaMP8f tau ~ 0.2-0.5s)
cfg.min_decay_tau_s = 0.05;  % Too fast = artifact
cfg.max_decay_tau_s = 8.0;   % Too slow = not returning to baseline

% =========================================================================
%  CORRELATION / DUPLICATE CRITERIA
% =========================================================================
cfg.max_pairwise_corr = 0.7;
cfg.corr_keep_higher_snr = true;
cfg.min_center_distance_px = 5;

% =========================================================================
%  ACTIVITY CRITERIA
% =========================================================================
cfg.min_event_count = 1;
cfg.max_snr = 150;  % Extreme artifacts only

% =========================================================================
%  BASELINE CRITERIA
% =========================================================================
cfg.require_baseline_return = true;
cfg.baseline_window_s = 5;
cfg.baseline_return_threshold = 0.5;

% =========================================================================
%  OUTPUT OPTIONS
% =========================================================================
cfg.generate_gallery = true;
cfg.generate_analysis_figures = true;
cfg.generate_rejection_examples = true;
cfg.gallery_segment_duration = 30;  % seconds per segment
cfg.skip_existing = false;          % Set false to rerun with new criteria
cfg.organize_subfolders = true;

end