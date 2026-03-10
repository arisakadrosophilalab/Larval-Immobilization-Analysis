function cfg = compare_conditions_config()
% COMPARE_CONDITIONS_CONFIG  Default configuration for COMPARE_NEURO_CONDITIONS_v3.
%
%   cfg = compare_conditions_config() returns a struct containing all
%   tunable parameters for the statistical comparison pipeline.
%
%   EXAMPLE:
%     cfg = compare_conditions_config();
%     cfg.conditions(1).name = 'Control';
%     cfg.conditions(1).refined_dir = '/path/to/refined/data';
%     cfg.conditions(2).name = 'Experimental';
%     cfg.conditions(2).refined_dir = '/path/to/refined/data';
%     cfg.output_dir = '/path/to/output';
%     COMPARE_NEURO_CONDITIONS_v3(cfg);
%
%   If conditions are not specified, a folder picker opens.
%
%   See also: COMPARE_NEURO_CONDITIONS_v3

% =========================================================================
%  DATA LOCATIONS
% =========================================================================
% Each condition needs:
%   .name         - Label (e.g., 'Control', 'Experimental')
%   .short_name   - Abbreviation (e.g., 'Ctrl', 'Exp')
%   .refined_dir  - Folder containing *_refined.mat files
%   .color        - RGB triplet for plotting
%   .color_light  - Lighter RGB for fill/shading
%   .raw_dir      - (optional) Folder with raw data for 'both' comparison mode
%
% If left empty, the pipeline will prompt for folder selection.
cfg.conditions = struct('name', {}, 'short_name', {}, 'refined_dir', {}, ...
    'color', {}, 'color_light', {}, 'raw_dir', {});

% Output directory for figures, tables, and reports.
cfg.output_dir = '';

% Path to neuroimaging motion data directory (from motion tracking).
% Set to '' to skip motion data integration.
cfg.neuro_motion_dir = '';

% =========================================================================
%  IMAGING PARAMETERS
% =========================================================================
cfg.frame_rate    = 33.33;   % Hz
cfg.pixel_size_um = 0.260;   % µm/pixel

% =========================================================================
%  STATISTICAL OPTIONS
% =========================================================================
cfg.stats = struct();
cfg.stats.alpha                = 0.05;     % Significance level
cfg.stats.bootstrap_iterations = 10000;    % For confidence intervals
cfg.stats.use_bonferroni       = true;     % Correct for multiple comparisons
cfg.stats.n_comparisons        = [];       % Counted dynamically if empty
cfg.stats.report_median        = true;     % Report median (for Mann-Whitney)
cfg.stats.compute_ci           = true;     % Compute 95% CI for effect sizes

% =========================================================================
%  COMPARISON MODE
% =========================================================================
% 'refined' - Only compare refined data (default, recommended for results)
% 'raw'     - Only compare raw data
% 'both'    - Compare both (recommended for methods papers)
cfg.comparison_mode = 'refined';

% =========================================================================
%  OUTPUT OPTIONS
% =========================================================================
cfg.output = struct();
cfg.output.generate_summary = true;  % Generate plain-language summary

% =========================================================================
%  FIGURE SETTINGS (Publication Quality)
% =========================================================================
cfg.fig = struct();
cfg.fig.font_name      = 'Arial';
cfg.fig.font_size_axis = 11;
cfg.fig.font_size_label = 13;
cfg.fig.font_size_title = 14;
cfg.fig.line_width     = 1.5;
cfg.fig.marker_size    = 40;
cfg.fig.dpi            = 300;
cfg.fig.save_formats   = {'png', 'pdf', 'svg'};

% =========================================================================
%  RANDOM SEED
% =========================================================================
cfg.rng_seed = 42;

end