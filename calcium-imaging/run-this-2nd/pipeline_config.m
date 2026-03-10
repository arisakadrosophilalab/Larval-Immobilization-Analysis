function config = pipeline_config()
% PIPELINE_CONFIG  Default configuration for COMPLETE_PIPELINE_RUN_ME_v4_5.
%
%   config = pipeline_config() returns a struct containing all tunable
%   parameters for the calcium imaging ROI detection pipeline. Modify
%   fields as needed before passing to the pipeline:
%
%   EXAMPLE:
%     config = pipeline_config();
%     config.pixelSize_um = 0.5;           % your pixel size in µm
%     config.frameRate = 30;                % your frame rate in Hz
%     config.conditions(1).name = 'Control';
%     config.conditions(1).tifDir = '/path/to/tif';
%     config.conditions(1).aviDir = '/path/to/avi';
%     COMPLETE_PIPELINE_RUN_ME_v4_5(config);
%
%   If conditions are not specified (or paths don't exist), the pipeline
%   will prompt for folder selection via a GUI picker.
%
%   See also: COMPLETE_PIPELINE_RUN_ME_v4_5

% =========================================================================
%  DATA LOCATIONS
% =========================================================================
% Define one or more conditions. Each condition needs:
%   .name    - Label (e.g., 'Control', 'Experimental')
%   .tifDir  - Folder containing raw TIF stacks
%   .aviDir  - Folder containing AVI files and shift CSVs from
%              analyze_motion_and_QC_neuro
%
% If left empty, the pipeline will prompt for folder selection.
config.conditions = struct('name', {}, 'tifDir', {}, 'aviDir', {});

% =========================================================================
%  ACQUISITION PARAMETERS
% =========================================================================
config.frameRate     = 1000/30;      % Hz (33.33 Hz for calcium recordings)
config.pixelSize_um  = 266/1024;     % µm/pixel (50x with 200 mm tube lens, 2x2 binning)
config.indicator     = 'GCaMP8f';    % Calcium indicator name (for provenance)

% Glob pattern for TIF file discovery. Only files matching this pattern
% in each TIF folder will be processed.
config.file_pattern  = 'S*.tif';     % e.g., '*.tif', 'S*.tif', 'fly*.tiff'

% =========================================================================
%  MOTION CORRECTION
% =========================================================================
config.motion = struct();
config.motion.chunk_size           = 500;   % Frames per chunk for motion correction
config.motion.max_reasonable_shift = 50;    % Warn if shifts exceed this (pixels)

% =========================================================================
%  ROI DETECTION (Drosophila L2 neuron sizes)
% =========================================================================
% Physical size calculations (pixel = 0.260 µm at default settings):
%   gSig = 3 px → radius = 0.78 µm → expected soma diameter ~1.6 µm
%   min_area = 12 px → ~0.8 µm² → min diameter ~1.0 µm
%   max_area = 120 px → ~8 µm² → max diameter ~3.2 µm
%   target_area = π×(3×1.2)² ≈ 41 px → ~2.8 µm² → diameter ~1.9 µm
%
% These are appropriate for Drosophila L2 neuron somata (2-3 µm diameter).
config.roi = struct();
config.roi.gSig              = 3;      % Expected neuron radius in pixels
config.roi.gSiz              = 8;      % ROI template size
config.roi.min_corr          = 0.3;    % Minimum Cn for detection
config.roi.min_pnr           = 2.0;    % Minimum PNR for detection
config.roi.min_snr           = 2.0;    % Minimum SNR for QC
config.roi.min_area_px       = 12;     % Minimum ROI area in pixels
config.roi.max_area_px       = 120;    % Maximum ROI area in pixels
config.roi.max_rois          = 300;    % Maximum ROIs per file
config.roi.neuropil_coef     = 0.7;    % Neuropil subtraction coefficient
config.roi.target_area_factor = 1.2;   % Growth factor for seed expansion
config.roi.local_percentile  = 30;     % Grow into top (100-X)% locally

% Global signal regression for widefield imaging
config.roi.apply_global_regression   = true;     % Remove shared neuropil signal
config.roi.global_regression_method  = 'linear'; % 'linear' or 'robust'

% =========================================================================
%  EVENT DETECTION (GCaMP8f optimized)
% =========================================================================
config.events = struct();
config.events.min_prominence  = 0.03;   % Minimum prominence (ΔF/F)
config.events.min_width_s     = 0.2;    % Minimum width (seconds)
config.events.max_width_s     = 10;     % Maximum width (seconds)
config.events.min_distance_s  = 0.3;    % Minimum inter-event distance (seconds)
config.events.min_rise_rate   = 0.01;   % Minimum rise rate (ΔF/F per second)

% =========================================================================
%  PROCESSING OPTIONS
% =========================================================================
config.skip_existing         = true;    % Skip already-processed files
config.save_shifts           = true;    % Save shift data
config.generate_figures      = true;    % Generate diagnostic figures
config.generate_trace_fig    = true;    % Generate calcium trace figures
config.generate_roi_map      = true;    % Generate clean ROI map
config.generate_roi_gallery  = true;    % Generate zoomed ROI thumbnails
config.save_diagnostic_figures = true;  % Save diagnostic figures
config.generate_pdfs         = true;    % Generate PDF versions of figures
config.force_regenerate_figs = false;   % Set true to regenerate ALL figures

% =========================================================================
%  FIGURE VISUALIZATION OPTIONS
% =========================================================================
config.fig = struct();
config.fig.roi_line_width     = 2.5;       % Outline line width
config.fig.roi_fill_alpha     = 0.3;       % Semi-transparent fill alpha
config.fig.roi_outline_color  = 'white';   % Outline color for contrast
config.fig.show_roi_numbers   = true;      % Label ROIs with numbers
config.fig.gallery_cols       = 6;         % Columns in ROI gallery
config.fig.gallery_rows       = 5;         % Rows per page in gallery

% =========================================================================
%  QC GATING
% =========================================================================
config.qc_gate_column       = 'ROI_ELIGIBLE';  % Column to gate on
config.require_qc_results   = true;            % Require QC results to proceed

% =========================================================================
%  RANDOM SEED
% =========================================================================
config.rng_seed = 42;

end