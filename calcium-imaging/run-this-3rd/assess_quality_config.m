function cfg = assess_quality_config()
% ASSESS_QUALITY_CONFIG  Default configuration for Assess_ROI_Quality.
%
%   cfg = assess_quality_config() returns a struct with all tunable
%   parameters for the post-pipeline quality assessment.
%
%   EXAMPLE:
%     cfg = assess_quality_config();
%     cfg.conditions(1).name = 'Control';
%     cfg.conditions(1).roi_dir = '/path/to/roi_analysis';
%     cfg.conditions(1).color = [0.2 0.4 0.8];
%     cfg.output_dir = '/path/to/output';
%     Assess_ROI_Quality(cfg);
%
%   If conditions are not specified, a folder picker opens.
%
%   See also: Assess_ROI_Quality

% =========================================================================
%  DATA LOCATIONS
% =========================================================================
% Each condition needs:
%   .name     - Label (e.g., 'Control', 'Experimental')
%   .roi_dir  - Folder containing *_rois.mat files from pipeline
%   .color    - RGB triplet for plotting
cfg.conditions = struct('name', {}, 'roi_dir', {}, 'color', {});

% Output directory for quality assessment results.
cfg.output_dir = '';

% =========================================================================
%  QUALITY THRESHOLDS (calibrated for pairwise synchrony)
% =========================================================================
cfg.thresholds = struct();
cfg.thresholds.max_dff          = 5;     % Artifact: extreme dF/F
cfg.thresholds.min_dff          = -2;    % Artifact: extreme negative dF/F
cfg.thresholds.severe_neuropil  = 0.50;  % Very high synchrony → unusable
cfg.thresholds.high_neuropil    = 0.30;  % Needs global regression
cfg.thresholds.mild_neuropil    = 0.15;  % Acceptable
cfg.thresholds.min_rois         = 5;     % Minimum ROIs for assessment

% =========================================================================
%  OUTPUT OPTIONS
% =========================================================================
cfg.options = struct();
cfg.options.generate_diagnostic_figures = true;
cfg.options.max_diagnostic_figures      = 10;
cfg.options.verbose                     = true;

end