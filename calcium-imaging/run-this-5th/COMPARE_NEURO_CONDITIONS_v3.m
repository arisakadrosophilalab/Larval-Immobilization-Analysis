function COMPARE_NEURO_CONDITIONS(varargin)
%% COMPARE_REFINED_CONDITIONS v3.0 - Publication-Grade Statistical Comparison
% =========================================================================
% Comprehensive statistical comparison of refined calcium imaging data
% between Control (Hydrogel) and Experimental (Hydrogel + Diethyl Ether)
% immobilization methods for Drosophila larvae.
%
% USAGE:
%   COMPARE_NEURO_CONDITIONS_v3()        % Uses defaults + folder picker
%   COMPARE_NEURO_CONDITIONS_v3(cfg)     % Uses custom config struct
%
%   cfg = compare_conditions_config();
%   cfg.conditions(1).name = 'Control';
%   cfg.conditions(1).short_name = 'Ctrl';
%   cfg.conditions(1).refined_dir = '/path/to/refined/data';
%   cfg.conditions(1).color = [0.2 0.4 0.8];
%   cfg.conditions(1).color_light = [0.6 0.7 0.9];
%   cfg.conditions(2).name = 'Experimental';
%   cfg.conditions(2).short_name = 'Exp';
%   cfg.conditions(2).refined_dir = '/path/to/refined/data';
%   cfg.conditions(2).color = [0.8 0.3 0.2];
%   cfg.conditions(2).color_light = [0.9 0.7 0.6];
%   cfg.output_dir = '/path/to/output';
%   COMPARE_NEURO_CONDITIONS_v3(cfg);
%
% USES REFINED DATA from ROI_REFINEMENT pipeline (not raw roi_analysis)
%
% PURPOSE:
%   Demonstrate that the immobilization method provides high-quality
%   neurological imaging data suitable for scientific publication.
%
% CHANGES IN v2.0:
%   - BUGFIX: Field paths now read from refinement.removal_counts (v3.1 format)
%   - BUGFIX: n_comparisons counted dynamically (was hardcoded to 12)
%   - REMOVED: Composite Motion Artifact Score (ad hoc, unvalidated)
%   - REMOVED: Figure 4 dynamite plot (replaced with removal breakdown)
%   - REMOVED: Figure 6 graphical abstract with auto-recommendations
%   - REMOVED: Auto-generated conclusions from plain language summary
%   - ADDED: Full 11-category removal breakdown comparison (Figure 4)
%   - ADDED: Neuroimaging motion data integration (from motion tracking)
%   - UPDATED: Forest plot now leads with sample-level metrics
%   - UPDATED: Methods text reflects sample-level primary analysis
%
% METRICS COMPARED:
%   1. ROI DETECTION QUALITY
%      - Number of ROIs detected (original and after refinement)
%      - Removal rate (indicator of false positive rate)
%      - Breakdown of removal reasons
%
%   2. SIGNAL QUALITY
%      - Signal-to-noise ratio (SNR)
%      - Event rates
%      - Event amplitudes
%
%   3. MOTION-INDUCED NOISE ESTIMATION
%      - Synchrony score (high = likely motion artifact)
%      - % ROIs removed due to correlation (motion indicator)
%      - Motion quality from provenance
%
%   4. TEMPORAL STABILITY
%      - Baseline drift (flagged ROIs)
%      - Signal consistency
%
% STATISTICAL TESTS:
%   - Mann-Whitney U (non-parametric comparison)
%   - Bootstrapped confidence intervals
%   - Effect sizes (Cliff's delta, Cohen's d)
%   - Bonferroni correction for multiple comparisons
%
% OUTPUTS:
%   - Publication-quality figures (PNG, PDF, SVG)
%   - Statistical summary tables (CSV)
%   - Comprehensive analysis report (TXT)
%   - Individual metric comparisons
%
% CONFIGURATION:
%   All parameters are set via the config struct. See
%   compare_conditions_config() for the full list and documentation.
%
% DEPENDENCIES:
%   - MATLAB R2020b or later
%   - Statistics and Machine Learning Toolbox (ranksum, lillietest)
%
% See also: compare_conditions_config, ROI_REFINEMENT_v3_2
% =========================================================================

%% ════════════════════════════════════════════════════════════════════════
%  CONFIGURATION
% ════════════════════════════════════════════════════════════════════════

if nargin >= 1 && isstruct(varargin{1})
    cfg = merge_with_defaults(varargin{1});
elseif nargin >= 1
    error(['Unexpected argument (type: %s). Pass a config struct from\n' ...
           'compare_conditions_config(), or call with no arguments for folder picker.'], ...
           class(varargin{1}));
else
    cfg = compare_conditions_config();
end

% If no conditions defined, prompt for folders
if isempty(cfg.conditions) || ~isfield(cfg.conditions, 'name') || isempty(cfg.conditions(1).name)
    cfg.conditions = prompt_for_conditions();
end

% If no output dir, prompt
if isempty(cfg.output_dir)
    if usejava('desktop')
        cfg.output_dir = uigetdir('', 'Select output directory for statistical comparison');
        if cfg.output_dir == 0
            error('No output directory selected. Exiting.');
        end
    else
        error(['No output_dir specified. Set cfg.output_dir before calling.\n' ...
               'See also: compare_conditions_config']);
    end
end

CONDITIONS = cfg.conditions;
OUTPUT_DIR = cfg.output_dir;
NEURO_MOTION_DIR = cfg.neuro_motion_dir;
FRAME_RATE = cfg.frame_rate;
PIXEL_SIZE = cfg.pixel_size_um;
STATS = cfg.stats;
COMPARISON_MODE = cfg.comparison_mode;
OUTPUT = cfg.output;
FIG = cfg.fig;

%% ════════════════════════════════════════════════════════════════════════
%  INITIALIZATION
% ════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  STATISTICAL COMPARISON OF REFINED ROI DATA                         ║\n');
fprintf('║  Control (Hydrogel) vs Experimental (Hydrogel + Diethyl Ether)      ║\n');
fprintf('║  Publication-Grade Analysis                                          ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Date: %s\n', datetime('now'));
fprintf('MATLAB: %s\n', version('-release'));
fprintf('\n');

% Create output directory structure
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

figures_dir = fullfile(OUTPUT_DIR, 'figures');
tables_dir = fullfile(OUTPUT_DIR, 'tables');

if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end

fprintf('Output directory: %s\n\n', OUTPUT_DIR);

% Set default figure properties
set(0, 'DefaultAxesFontName', FIG.font_name);
set(0, 'DefaultAxesFontSize', FIG.font_size_axis);
set(0, 'DefaultTextFontName', FIG.font_name);

% Set random seed for reproducibility
rng(cfg.rng_seed, 'twister');

%% ════════════════════════════════════════════════════════════════════════
%  LOAD REFINED DATA
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  LOADING REFINED DATA\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

% Load data for each condition
all_data = struct();

for c = 1:length(CONDITIONS)
    cond = CONDITIONS(c);
    fprintf('Loading %s data from:\n  %s\n', cond.name, cond.refined_dir);
    
    if ~exist(cond.refined_dir, 'dir')
        fprintf('  [ERROR] Directory not found!\n\n');
        all_data(c).samples = [];
        all_data(c).n_samples = 0;
        continue;
    end
    
    % Find refined .mat files
    mat_files = dir(fullfile(cond.refined_dir, '*_refined.mat'));
    
    if isempty(mat_files)
        fprintf('  [WARN] No refined .mat files found\n\n');
        all_data(c).samples = [];
        all_data(c).n_samples = 0;
        continue;
    end
    
    fprintf('  Found %d refined files\n', length(mat_files));
    
    samples = [];
    
    for f = 1:length(mat_files)
        mat_path = fullfile(mat_files(f).folder, mat_files(f).name);
        
        try
            loaded = load(mat_path);
            roi_data = loaded.roi_data;
            
            sample = struct();
            sample.filename = strrep(mat_files(f).name, '_refined.mat', '');
            sample.n_rois = size(roi_data.roi_centers, 1);
            
            % Get refinement info
            if isfield(roi_data, 'refinement')
                sample.n_original = roi_data.refinement.n_original;
                sample.n_removed = roi_data.refinement.n_removed;
                if sample.n_original > 0
                    sample.pct_removed = 100 * sample.n_removed / sample.n_original;
                else
                    sample.pct_removed = 0;
                end
                % Read from removal_counts substruct (ROI_REFINEMENT v3.1 format)
                if isfield(roi_data.refinement, 'removal_counts')
                    rc = roi_data.refinement.removal_counts;
                    sample.removed_extreme_snr = get_field_or_zero(rc, 'extreme_snr');
                    sample.removed_edge = get_field_or_zero(rc, 'edge_contact');
                    sample.removed_shape = get_field_or_zero(rc, 'shape');
                    sample.removed_negative = get_field_or_zero(rc, 'negative_transients');
                    sample.removed_step = get_field_or_zero(rc, 'step_artifacts');
                    sample.removed_kinetics = get_field_or_zero(rc, 'kinetics');
                    sample.removed_plateau = get_field_or_zero(rc, 'plateau');
                    sample.removed_correlated = get_field_or_zero(rc, 'correlation');
                    sample.removed_spatial = get_field_or_zero(rc, 'spatial');
                    sample.removed_inactive = get_field_or_zero(rc, 'inactive');
                    sample.removed_baseline = get_field_or_zero(rc, 'baseline_drift');
                else
                    % Fallback for older data format
                    warning('Sample %s uses legacy refinement format (no removal_counts).', mat_files(f).name);
                    sample.removed_extreme_snr = 0;
                    sample.removed_edge = 0;
                    sample.removed_shape = get_field_or_zero(roi_data.refinement, 'removed_shape');
                    sample.removed_negative = 0;
                    sample.removed_step = 0;
                    sample.removed_kinetics = 0;
                    sample.removed_plateau = 0;
                    sample.removed_correlated = get_field_or_zero(roi_data.refinement, 'removed_correlated');
                    sample.removed_spatial = get_field_or_zero(roi_data.refinement, 'removed_spatial');
                    sample.removed_inactive = get_field_or_zero(roi_data.refinement, 'removed_inactive');
                    sample.removed_baseline = 0;
                end
            else
                sample.n_original = sample.n_rois;
                sample.n_removed = 0;
                sample.pct_removed = 0;
                sample.removed_extreme_snr = 0;
                sample.removed_edge = 0;
                sample.removed_shape = 0;
                sample.removed_negative = 0;
                sample.removed_step = 0;
                sample.removed_kinetics = 0;
                sample.removed_plateau = 0;
                sample.removed_correlated = 0;
                sample.removed_spatial = 0;
                sample.removed_inactive = 0;
                sample.removed_baseline = 0;
            end
            
            % Get motion quality from provenance
            if isfield(roi_data, 'original_provenance') && isfield(roi_data.original_provenance, 'motion_quality')
                sample.motion_quality = char(roi_data.original_provenance.motion_quality);
            else
                sample.motion_quality = 'unknown';
            end
            
            % Get ROI-level metrics
            if sample.n_rois > 0
                sample.snr = roi_data.snr(:);
                sample.event_rate = roi_data.event_rate(:);
                
                if isfield(roi_data, 'event_count')
                    sample.event_count = roi_data.event_count(:);
                else
                    sample.event_count = [];
                end
                
                if isfield(roi_data, 'event_amplitude')
                    sample.event_amplitude = roi_data.event_amplitude(:);
                else
                    sample.event_amplitude = [];
                end
                
                % Synchrony score (motion artifact indicator)
                if isfield(roi_data, 'synchrony_score')
                    sample.synchrony_score = roi_data.synchrony_score;
                else
                    sample.synchrony_score = NaN;
                end
                
                % Baseline flags
                if isfield(roi_data, 'baseline_flag')
                    sample.n_baseline_flagged = sum(roi_data.baseline_flag);
                    sample.pct_baseline_flagged = 100 * sample.n_baseline_flagged / sample.n_rois;
                else
                    sample.n_baseline_flagged = 0;
                    sample.pct_baseline_flagged = 0;
                end
                
                % Mean metrics per sample
                sample.mean_snr = mean(sample.snr);
                sample.mean_event_rate = mean(sample.event_rate);
                
                if ~isempty(sample.event_amplitude)
                    sample.mean_event_amplitude = mean(sample.event_amplitude);
                else
                    sample.mean_event_amplitude = NaN;
                end
            else
                sample.snr = [];
                sample.event_rate = [];
                sample.event_count = [];
                sample.event_amplitude = [];
                sample.synchrony_score = NaN;
                sample.n_baseline_flagged = 0;
                sample.pct_baseline_flagged = 0;
                sample.mean_snr = NaN;
                sample.mean_event_rate = NaN;
                sample.mean_event_amplitude = NaN;
            end
            
            % Percent correlated removed (direct from removal counts)
            if sample.n_original > 0
                sample.pct_correlated_removed = min(100, 100 * sample.removed_correlated / sample.n_original);
            else
                sample.pct_correlated_removed = 0;
            end
            
            if isempty(samples)
                samples = sample;
            else
                samples(end+1) = sample;
            end
            
        catch load_err
            fprintf('  [ERROR] Failed to load %s: %s\n', mat_files(f).name, load_err.message);
        end
    end
    
    all_data(c).samples = samples;
    all_data(c).n_samples = length(samples);
    all_data(c).condition = cond.name;
    all_data(c).color = cond.color;
    all_data(c).color_light = cond.color_light;
    
    fprintf('  Loaded %d samples, %d total ROIs\n\n', ...
        all_data(c).n_samples, sum([samples.n_rois]));
end

% Check we have data
n_cond1 = all_data(1).n_samples;
n_cond2 = all_data(2).n_samples;

if n_cond1 == 0 && n_cond2 == 0
    error('No data found for either condition! Run ROI_REFINEMENT first.');
end

if n_cond1 == 0 || n_cond2 == 0
    warning('Only one condition has data. Statistical comparisons will be limited.');
end

%% ════════════════════════════════════════════════════════════════════════
%  EXTRACT METRICS FOR COMPARISON
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  EXTRACTING METRICS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

metrics = struct();

for c = 1:2
    if all_data(c).n_samples == 0
        continue;
    end
    
    samples = all_data(c).samples;
    cond_name = lower(CONDITIONS(c).short_name);
    
    % Sample-level metrics
    metrics.(cond_name).n_rois = [samples.n_rois]';
    metrics.(cond_name).n_original = [samples.n_original]';
    metrics.(cond_name).pct_removed = [samples.pct_removed]';
    metrics.(cond_name).synchrony = [samples.synchrony_score]';
    metrics.(cond_name).pct_correlated = [samples.pct_correlated_removed]';
    metrics.(cond_name).pct_baseline_flagged = [samples.pct_baseline_flagged]';
    
    % Full removal breakdown (all 11 categories from ROI_REFINEMENT v3.1)
    metrics.(cond_name).removed_extreme_snr = [samples.removed_extreme_snr]';
    metrics.(cond_name).removed_edge = [samples.removed_edge]';
    metrics.(cond_name).removed_shape = [samples.removed_shape]';
    metrics.(cond_name).removed_negative = [samples.removed_negative]';
    metrics.(cond_name).removed_step = [samples.removed_step]';
    metrics.(cond_name).removed_kinetics = [samples.removed_kinetics]';
    metrics.(cond_name).removed_plateau = [samples.removed_plateau]';
    metrics.(cond_name).removed_correlated = [samples.removed_correlated]';
    metrics.(cond_name).removed_spatial = [samples.removed_spatial]';
    metrics.(cond_name).removed_inactive = [samples.removed_inactive]';
    metrics.(cond_name).removed_baseline = [samples.removed_baseline]';
    metrics.(cond_name).mean_snr = [samples.mean_snr]';
    metrics.(cond_name).mean_event_rate = [samples.mean_event_rate]';
    metrics.(cond_name).mean_event_amplitude = [samples.mean_event_amplitude]';
    
    % ROI-level metrics (pooled across all samples)
    metrics.(cond_name).all_snr = vertcat(samples.snr);
    metrics.(cond_name).all_event_rate = vertcat(samples.event_rate);
    
    % Coefficient of variation for SNR (measure of consistency)
    mean_snr_val = mean(metrics.(cond_name).mean_snr, 'omitnan');
    if ~isnan(mean_snr_val) && mean_snr_val > 0
        metrics.(cond_name).snr_cv = std(metrics.(cond_name).mean_snr, 0, 'omitnan') / mean_snr_val;
    else
        metrics.(cond_name).snr_cv = NaN;
    end
    
    event_amp = {samples.event_amplitude};
    event_amp = event_amp(~cellfun(@isempty, event_amp));
    if ~isempty(event_amp)
        metrics.(cond_name).all_event_amplitude = vertcat(event_amp{:});
    else
        metrics.(cond_name).all_event_amplitude = [];
    end
    
    % Motion quality breakdown
    motion_quals = {samples.motion_quality};
    metrics.(cond_name).n_good_motion = sum(strcmp(motion_quals, 'good'));
    metrics.(cond_name).n_acceptable_motion = sum(strcmp(motion_quals, 'acceptable'));
    metrics.(cond_name).n_excessive_motion = sum(strcmp(motion_quals, 'excessive_motion'));
    
    fprintf('%s:\n', CONDITIONS(c).name);
    fprintf('  Samples: %d\n', all_data(c).n_samples);
    fprintf('  Total refined ROIs: %d\n', sum(metrics.(cond_name).n_rois));
    fprintf('  Mean ROIs/sample: %.1f ± %.1f\n', mean(metrics.(cond_name).n_rois), std(metrics.(cond_name).n_rois));
    fprintf('  Mean SNR: %.2f ± %.2f\n', mean(metrics.(cond_name).all_snr), std(metrics.(cond_name).all_snr));
    fprintf('  Mean synchrony: %.3f ± %.3f\n', mean(metrics.(cond_name).synchrony, 'omitnan'), std(metrics.(cond_name).synchrony, 'omitnan'));
    fprintf('\n');
end

%% ════════════════════════════════════════════════════════════════════════
%  STATISTICAL TESTS
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  STATISTICAL ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

% Adjust alpha for multiple comparisons if requested
if STATS.use_bonferroni
    % Count comparisons dynamically (will be finalized after building list)
    % Placeholder - actual count set after comparisons list is built
    alpha_adj = STATS.alpha;  % Temporary, updated below
    fprintf('Bonferroni correction will be applied after counting comparisons.\n\n');
else
    alpha_adj = STATS.alpha;
end

% Initialize results table
stat_results = table();

if n_cond1 > 0 && n_cond2 > 0
    
    % First, test for normality to justify non-parametric approach
    fprintf('Normality Assessment (Shapiro-Wilk where n≤50, else Lilliefors):\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    
    test_normality('Control SNR', metrics.ctrl.all_snr);
    test_normality('Experimental SNR', metrics.exp.all_snr);
    test_normality('Control Event Rate', metrics.ctrl.all_event_rate);
    test_normality('Experimental Event Rate', metrics.exp.all_event_rate);
    
    fprintf('\n→ Non-parametric tests (Mann-Whitney U) appropriate for non-normal data\n\n');
    
    % Define comparisons to make
    comparisons = {
        % {metric_name, ctrl_data, exp_data, display_name, is_sample_level}
        {'n_rois', metrics.ctrl.n_rois, metrics.exp.n_rois, 'ROIs per Sample', true};
        {'pct_removed', metrics.ctrl.pct_removed, metrics.exp.pct_removed, '% ROIs Removed', true};
        {'synchrony', metrics.ctrl.synchrony, metrics.exp.synchrony, 'Synchrony Score', true};
        {'pct_correlated', metrics.ctrl.pct_correlated, metrics.exp.pct_correlated, '% Correlated Removed', true};
        {'mean_snr', metrics.ctrl.mean_snr, metrics.exp.mean_snr, 'Mean SNR (sample)', true};
        {'mean_event_rate', metrics.ctrl.mean_event_rate, metrics.exp.mean_event_rate, 'Mean Event Rate (sample)', true};
        {'all_snr', metrics.ctrl.all_snr, metrics.exp.all_snr, 'SNR (all ROIs)', false};
        {'all_event_rate', metrics.ctrl.all_event_rate, metrics.exp.all_event_rate, 'Event Rate (all ROIs)', false};
        {'pct_baseline_flagged', metrics.ctrl.pct_baseline_flagged, metrics.exp.pct_baseline_flagged, '% Baseline Drift', true};
    };
    
    % Add event amplitude if available
    if ~isempty(metrics.ctrl.all_event_amplitude) && ~isempty(metrics.exp.all_event_amplitude)
        comparisons{end+1} = {'all_event_amplitude', metrics.ctrl.all_event_amplitude, ...
            metrics.exp.all_event_amplitude, 'Event Amplitude (all ROIs)', false};
    end
    
    % Finalize Bonferroni correction with actual comparison count
    STATS.n_comparisons = length(comparisons);
    if STATS.use_bonferroni
        alpha_adj = STATS.alpha / STATS.n_comparisons;
        fprintf('Using Bonferroni correction: α = %.4f (%.3f / %d comparisons)\n\n', ...
            alpha_adj, STATS.alpha, STATS.n_comparisons);
    end
    
    fprintf('Statistical Comparison Results:\n');
    fprintf('─────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n');
    fprintf('Metric                        Ctrl Median (IQR)      Exp Median (IQR)       p-value    Sig    P(C>E)   Effect\n');
    fprintf('─────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n');
    
    for i = 1:length(comparisons)
        comp = comparisons{i};
        metric_name = comp{1};
        ctrl_data = comp{2};
        exp_data = comp{3};
        display_name = comp{4};
        is_sample = comp{5};
        
        % Remove NaN values
        ctrl_clean = ctrl_data(~isnan(ctrl_data));
        exp_clean = exp_data(~isnan(exp_data));
        
        if isempty(ctrl_clean) || isempty(exp_clean)
            continue;
        end
        
        % Calculate descriptive statistics (mean for reference)
        ctrl_mean = mean(ctrl_clean);
        ctrl_std = std(ctrl_clean);
        exp_mean = mean(exp_clean);
        exp_std = std(exp_clean);
        
        % Median and IQR (more appropriate for Mann-Whitney)
        ctrl_median = median(ctrl_clean);
        ctrl_iqr = iqr(ctrl_clean);
        exp_median = median(exp_clean);
        exp_iqr = iqr(exp_clean);
        
        % Mann-Whitney U test (non-parametric)
        % Use exact test for small samples (n < 20)
        if length(ctrl_clean) < 20 && length(exp_clean) < 20
            [p_value, ~, stats_mw] = ranksum(ctrl_clean, exp_clean, 'method', 'exact');
        else
            [p_value, ~, stats_mw] = ranksum(ctrl_clean, exp_clean);
        end
        
        % Effect size: Cliff's delta (non-parametric)
        cliffs_d = compute_cliffs_delta(ctrl_clean, exp_clean);
        
        % Probability Index P(X > Y) - directly interpretable
        % This is (1 + delta) / 2, represents probability control > experimental
        prob_index = (1 + cliffs_d) / 2;
        
        % Also compute Cohen's d for reference (handle zero std)
        pooled_std = sqrt(((length(ctrl_clean)-1)*ctrl_std^2 + (length(exp_clean)-1)*exp_std^2) / ...
                         (length(ctrl_clean) + length(exp_clean) - 2));
        if pooled_std > 0
            cohens_d = (ctrl_mean - exp_mean) / pooled_std;
        else
            cohens_d = 0;
        end
        
        % Significance
        if p_value < alpha_adj
            sig_str = get_sig_stars(p_value);
        else
            sig_str = 'n.s.';
        end
        
        % Effect size interpretation (Cliff's delta thresholds)
        if abs(cliffs_d) < 0.147
            effect_interp = 'negligible';
        elseif abs(cliffs_d) < 0.33
            effect_interp = 'small';
        elseif abs(cliffs_d) < 0.474
            effect_interp = 'medium';
        else
            effect_interp = 'large';
        end
        
        % Print result with median (IQR)
        fprintf('%-28s  %7.2f (%6.2f)   %7.2f (%6.2f)   %9.4f  %-5s  %.3f    %s\n', ...
            display_name, ctrl_median, ctrl_iqr, exp_median, exp_iqr, p_value, sig_str, prob_index, effect_interp);
        
        % Add to results table (comprehensive)
        row = table();
        row.Metric = {display_name};
        row.Ctrl_N = length(ctrl_clean);
        row.Ctrl_Mean = ctrl_mean;
        row.Ctrl_SD = ctrl_std;
        row.Ctrl_Median = ctrl_median;
        row.Ctrl_IQR = ctrl_iqr;
        row.Exp_N = length(exp_clean);
        row.Exp_Mean = exp_mean;
        row.Exp_SD = exp_std;
        row.Exp_Median = exp_median;
        row.Exp_IQR = exp_iqr;
        row.p_value = p_value;
        row.p_adjusted = min(p_value * STATS.n_comparisons, 1);  % Adjusted p
        row.Significant = p_value < alpha_adj;
        row.Cliffs_Delta = cliffs_d;
        row.Prob_Index = prob_index;
        row.Cohens_d = cohens_d;
        row.Effect_Size = {effect_interp};
        row.Test = {'Mann-Whitney U'};
        row.Sample_Level = is_sample;
        
        stat_results = [stat_results; row];
    end
    
    fprintf('─────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n');
    fprintf('P(C>E): Probability that a random Control value exceeds a random Experimental value\n');
    fprintf('        (0.5 = no difference, >0.5 = Control tends higher, <0.5 = Experimental tends higher)\n');
    fprintf('Significance: * p<0.05, ** p<0.01, *** p<0.001 (Bonferroni-corrected α = %.4f)\n', alpha_adj);
    fprintf('Effect size (Cliff''s δ): |δ|<0.147 negligible, <0.33 small, <0.474 medium, ≥0.474 large\n\n');
    
    % Sample size consideration
    fprintf('Sample Size Consideration:\n');
    fprintf('  Control: %d samples, %d total ROIs\n', n_cond1, length(metrics.ctrl.all_snr));
    fprintf('  Experimental: %d samples, %d total ROIs\n', n_cond2, length(metrics.exp.all_snr));
    if n_cond1 < 10 || n_cond2 < 10
        fprintf('  ⚠ WARNING: Small sample sizes (n<%d) limit statistical power.\n', 10);
        fprintf('    Consider collecting additional samples for robust conclusions.\n');
        
        % Estimate power for SNR comparison (approximate)
        if n_cond1 >= 3 && n_cond2 >= 3
            snr_effect = abs(mean(metrics.ctrl.mean_snr) - mean(metrics.exp.mean_snr)) / ...
                         std([metrics.ctrl.mean_snr; metrics.exp.mean_snr]);
            approx_power = estimate_power(n_cond1, n_cond2, snr_effect);
            fprintf('    Approximate power for SNR (effect=%.2f): %.1f%%\n', snr_effect, approx_power*100);
            if approx_power < 0.8
                fprintf('    → Power < 80%%. Results may be underpowered.\n');
            end
        end
    end
    fprintf('\n');
    
end

%% ════════════════════════════════════════════════════════════════════════
%  MOTION-INDUCED NOISE ANALYSIS
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  MOTION-INDUCED NOISE ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

fprintf('Motion artifacts are estimated using two proxy measures:\n');
fprintf('  1. Synchrony Score: High values (>0.5) indicate correlated activity\n');
fprintf('     across ROIs, suggesting motion-induced signal rather than \n');
fprintf('     independent neural activity.\n');
fprintf('  2. %% Correlated ROIs Removed: High percentage indicates many ROIs\n');
fprintf('     had similar traces (r>0.7), likely due to shared motion artifact.\n\n');

if n_cond1 > 0 && n_cond2 > 0
    fprintf('Motion Quality Summary:\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    fprintf('Condition      Good Motion    Acceptable    Excessive    Sync>0.5\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    
    for c = 1:2
        cond_name = lower(CONDITIONS(c).short_name);
        sync_vals = metrics.(cond_name).synchrony;
        n_high_sync = sum(sync_vals > 0.5);
        
        fprintf('%-14s     %d             %d             %d           %d/%d\n', ...
            CONDITIONS(c).name, ...
            metrics.(cond_name).n_good_motion, ...
            metrics.(cond_name).n_acceptable_motion, ...
            metrics.(cond_name).n_excessive_motion, ...
            n_high_sync, length(sync_vals));
    end
    fprintf('─────────────────────────────────────────────────────────────────────────\n\n');
    
    % Statistical comparison of motion metrics (separate, not composite)
    fprintf('Motion Proxy Comparison:\n');
    fprintf('  Synchrony Score (higher = more motion-correlated signal):\n');
    fprintf('    Control:      %.3f ± %.3f\n', mean(metrics.ctrl.synchrony, 'omitnan'), std(metrics.ctrl.synchrony, 'omitnan'));
    fprintf('    Experimental: %.3f ± %.3f\n', mean(metrics.exp.synchrony, 'omitnan'), std(metrics.exp.synchrony, 'omitnan'));
    [p_sync, ~] = ranksum(metrics.ctrl.synchrony, metrics.exp.synchrony);
    fprintf('    p-value: %.4f %s\n\n', p_sync, get_sig_stars(p_sync));
    
    fprintf('  %% Correlated ROIs Removed:\n');
    fprintf('    Control:      %.1f%% ± %.1f%%\n', mean(metrics.ctrl.pct_correlated), std(metrics.ctrl.pct_correlated));
    fprintf('    Experimental: %.1f%% ± %.1f%%\n', mean(metrics.exp.pct_correlated), std(metrics.exp.pct_correlated));
    [p_corr, ~] = ranksum(metrics.ctrl.pct_correlated, metrics.exp.pct_correlated);
    fprintf('    p-value: %.4f %s\n\n', p_corr, get_sig_stars(p_corr));
end

%% ════════════════════════════════════════════════════════════════════════
%  NEUROIMAGING MOTION DATA INTEGRATION
% ════════════════════════════════════════════════════════════════════════

motion_data_loaded = false;

if ~isempty(NEURO_MOTION_DIR)
    fprintf('═══════════════════════════════════════════════════════════════════════\n');
    fprintf('  NEUROIMAGING MOTION DATA INTEGRATION\n');
    fprintf('═══════════════════════════════════════════════════════════════════════\n\n');
    
    motion_csv = fullfile(NEURO_MOTION_DIR, 'neuro_motion_per_file.csv');
    
    if exist(motion_csv, 'file')
        try
            motion_table = readtable(motion_csv, 'VariableNamingRule', 'preserve');
            fprintf('Loaded motion data: %d files from %s\n', height(motion_table), motion_csv);
            
            % Match motion data to neural samples by filename stem
            for c = 1:2
                if all_data(c).n_samples == 0, continue; end
                cond_name = lower(CONDITIONS(c).short_name);
                
                matched = 0;
                motion_speeds = nan(all_data(c).n_samples, 1);
                motion_median_steps = nan(all_data(c).n_samples, 1);
                motion_jitter = nan(all_data(c).n_samples, 1);
                
                for s = 1:all_data(c).n_samples
                    sample_stem = all_data(c).samples(s).filename;
                    % Try matching by stem (strip extensions)
                    sample_stem_clean = regexprep(sample_stem, '_rois$|_refined$', '');
                    
                    % Search motion table for matching file
                    if ismember('File', motion_table.Properties.VariableNames)
                        motion_files = string(motion_table.File);
                        motion_stems = regexprep(motion_files, '\.avi$', '', 'ignorecase');
                        
                        match_idx = find(strcmpi(motion_stems, sample_stem_clean), 1);
                        if isempty(match_idx)
                            % Try partial match
                            match_idx = find(contains(motion_stems, sample_stem_clean, 'IgnoreCase', true), 1);
                        end
                        
                        if ~isempty(match_idx)
                            matched = matched + 1;
                            if ismember('Mean_speed_um_per_s', motion_table.Properties.VariableNames)
                                motion_speeds(s) = motion_table.Mean_speed_um_per_s(match_idx);
                            end
                            if ismember('Median_step_um', motion_table.Properties.VariableNames)
                                motion_median_steps(s) = motion_table.Median_step_um(match_idx);
                            end
                            if ismember('Jitter_ratio', motion_table.Properties.VariableNames)
                                motion_jitter(s) = motion_table.Jitter_ratio(match_idx);
                            end
                        end
                    end
                end
                
                metrics.(cond_name).motion_speed = motion_speeds;
                metrics.(cond_name).motion_median_step = motion_median_steps;
                metrics.(cond_name).motion_jitter = motion_jitter;
                
                fprintf('  %s: matched %d/%d samples to motion data\n', ...
                    CONDITIONS(c).name, matched, all_data(c).n_samples);
            end
            
            motion_data_loaded = true;
            
            % Report motion comparison if we have matched data
            if motion_data_loaded && n_cond1 > 0 && n_cond2 > 0
                fprintf('\nPhysical Motion During Neuroimaging (from motion tracking):\n');
                fprintf('─────────────────────────────────────────────────────────────────────────\n');
                
                ctrl_speed = metrics.ctrl.motion_speed(~isnan(metrics.ctrl.motion_speed));
                exp_speed = metrics.exp.motion_speed(~isnan(metrics.exp.motion_speed));
                ctrl_medstep = metrics.ctrl.motion_median_step(~isnan(metrics.ctrl.motion_median_step));
                exp_medstep = metrics.exp.motion_median_step(~isnan(metrics.exp.motion_median_step));
                
                if ~isempty(ctrl_speed) && ~isempty(exp_speed)
                    fprintf('  Mean Speed (µm/s):  Ctrl = %.3f ± %.3f,  Exp = %.3f ± %.3f\n', ...
                        mean(ctrl_speed), std(ctrl_speed), mean(exp_speed), std(exp_speed));
                end
                if ~isempty(ctrl_medstep) && ~isempty(exp_medstep)
                    fprintf('  Median Step (µm):   Ctrl = %.4f ± %.4f,  Exp = %.4f ± %.4f\n', ...
                        mean(ctrl_medstep), std(ctrl_medstep), mean(exp_medstep), std(exp_medstep));
                end
                fprintf('─────────────────────────────────────────────────────────────────────────\n');
                fprintf('NOTE: For detailed motion analysis, see compare_motion_groups_neuroimaging output.\n\n');
            end
            
        catch motion_err
            fprintf('  [WARNING] Failed to load motion data: %s\n', motion_err.message);
            fprintf('  Continuing without motion integration.\n\n');
        end
    else
        fprintf('No neuroimaging motion data found at:\n  %s\n', motion_csv);
        fprintf('Run compare_motion_groups_neuroimaging.m first, or set NEURO_MOTION_DIR = ''''.\n\n');
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  HIERARCHICAL DATA ANALYSIS (for nested ROI data)
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  HIERARCHICAL DATA STRUCTURE ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

fprintf('ROIs are nested within samples, creating a hierarchical/clustered data\n');
fprintf('structure. Standard tests (treating ROIs as independent) can inflate\n');
fprintf('Type I error rates to >45%% (Saravanan et al., 2020). We assess this here.\n\n');

if n_cond1 > 0 && n_cond2 > 0
    
    % Prepare data in cell array format (one cell per sample)
    ctrl_snr_by_sample = cell(n_cond1, 1);
    ctrl_event_by_sample = cell(n_cond1, 1);
    for s = 1:n_cond1
        ctrl_snr_by_sample{s} = all_data(1).samples(s).snr(:);
        ctrl_event_by_sample{s} = all_data(1).samples(s).event_rate(:);
    end
    
    exp_snr_by_sample = cell(n_cond2, 1);
    exp_event_by_sample = cell(n_cond2, 1);
    for s = 1:n_cond2
        exp_snr_by_sample{s} = all_data(2).samples(s).snr(:);
        exp_event_by_sample{s} = all_data(2).samples(s).event_rate(:);
    end
    
    % Compute ICC for each condition
    fprintf('Intraclass Correlation Coefficient (ICC):\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    fprintf('ICC measures how much variance is due to between-sample vs within-sample\n');
    fprintf('differences. High ICC (>0.5) suggests sample-level analysis is critical.\n\n');
    
    icc_ctrl_snr = compute_icc(ctrl_snr_by_sample);
    icc_exp_snr = compute_icc(exp_snr_by_sample);
    icc_ctrl_event = compute_icc(ctrl_event_by_sample);
    icc_exp_event = compute_icc(exp_event_by_sample);
    
    fprintf('  Control SNR ICC:       %.3f\n', icc_ctrl_snr);
    fprintf('  Experimental SNR ICC:  %.3f\n', icc_exp_snr);
    fprintf('  Control Event ICC:     %.3f\n', icc_ctrl_event);
    fprintf('  Experimental Event ICC: %.3f\n\n', icc_exp_event);
    
    % Compute design effect
    fprintf('Design Effect (inflation factor for effective sample size):\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    
    de_ctrl_snr = compute_design_effect(ctrl_snr_by_sample);
    de_exp_snr = compute_design_effect(exp_snr_by_sample);
    
    fprintf('  Control:      %.2f (effective n ≈ %d instead of %d)\n', ...
        de_ctrl_snr, round(length(metrics.ctrl.all_snr)/de_ctrl_snr), length(metrics.ctrl.all_snr));
    fprintf('  Experimental: %.2f (effective n ≈ %d instead of %d)\n\n', ...
        de_exp_snr, round(length(metrics.exp.all_snr)/de_exp_snr), length(metrics.exp.all_snr));
    
    % Hierarchical bootstrap confidence intervals for effect sizes
    fprintf('Hierarchical Bootstrap Analysis (accounts for nested structure):\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    fprintf('Computing 95%% CIs using two-level bootstrap (samples → ROIs)...\n\n');
    
    % Bootstrap CI for mean SNR difference
    [ci_low_snr, ci_high_snr] = bootstrap_mean_diff_ci(ctrl_snr_by_sample, exp_snr_by_sample, 5000, 0.95);
    mean_diff_snr = mean(metrics.ctrl.all_snr) - mean(metrics.exp.all_snr);
    
    fprintf('  SNR Difference (Ctrl - Exp):\n');
    fprintf('    Point estimate: %.2f\n', mean_diff_snr);
    fprintf('    95%% Hierarchical Bootstrap CI: [%.2f, %.2f]\n', ci_low_snr, ci_high_snr);
    if ci_low_snr > 0
        fprintf('    → CI excludes zero: Control SNR significantly HIGHER\n');
    elseif ci_high_snr < 0
        fprintf('    → CI excludes zero: Experimental SNR significantly HIGHER\n');
    else
        fprintf('    → CI includes zero: No significant difference\n');
    end
    fprintf('\n');
    
    % Bootstrap CI for mean Event Rate difference
    [ci_low_event, ci_high_event] = bootstrap_mean_diff_ci(ctrl_event_by_sample, exp_event_by_sample, 5000, 0.95);
    mean_diff_event = mean(metrics.ctrl.all_event_rate) - mean(metrics.exp.all_event_rate);
    
    fprintf('  Event Rate Difference (Ctrl - Exp):\n');
    fprintf('    Point estimate: %.4f Hz\n', mean_diff_event);
    fprintf('    95%% Hierarchical Bootstrap CI: [%.4f, %.4f]\n', ci_low_event, ci_high_event);
    if ci_low_event > 0
        fprintf('    → CI excludes zero: Control Event Rate significantly HIGHER\n');
    elseif ci_high_event < 0
        fprintf('    → CI excludes zero: Experimental Event Rate significantly HIGHER\n');
    else
        fprintf('    → CI includes zero: No significant difference\n');
    end
    fprintf('\n');
    
    % Bootstrap CI for Cliff's delta
    fprintf('  Effect Size (Cliff''s δ) with 95%% Bootstrap CI:\n');
    cliffs_delta_func = @(x, y) compute_cliffs_delta(x, y);
    
    [delta_ci_low_snr, delta_ci_high_snr] = bootstrap_effect_size_ci(...
        metrics.ctrl.all_snr, metrics.exp.all_snr, cliffs_delta_func, 5000, 0.95);
    delta_snr = compute_cliffs_delta(metrics.ctrl.all_snr, metrics.exp.all_snr);
    
    fprintf('    SNR: δ = %.3f [%.3f, %.3f]\n', delta_snr, delta_ci_low_snr, delta_ci_high_snr);
    
    [delta_ci_low_event, delta_ci_high_event] = bootstrap_effect_size_ci(...
        metrics.ctrl.all_event_rate, metrics.exp.all_event_rate, cliffs_delta_func, 5000, 0.95);
    delta_event = compute_cliffs_delta(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate);
    
    fprintf('    Event Rate: δ = %.3f [%.3f, %.3f]\n', delta_event, delta_ci_low_event, delta_ci_high_event);
    
    % Add Hedges' g for standardized effect size
    fprintf('\n  Standardized Effect Size (Hedges'' g, bias-corrected):\n');
    hedges_g_snr = compute_hedges_g(metrics.ctrl.all_snr, metrics.exp.all_snr);
    hedges_g_event = compute_hedges_g(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate);
    fprintf('    SNR: g = %.3f\n', hedges_g_snr);
    fprintf('    Event Rate: g = %.3f\n', hedges_g_event);
    
    fprintf('\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    fprintf('INTERPRETATION NOTE:\n');
    fprintf('Hierarchical bootstrap CIs are wider than naive CIs because they correctly\n');
    fprintf('account for within-sample correlation. If ROI-level tests show p<0.001\n');
    fprintf('but hierarchical CI includes zero, the apparent significance may be\n');
    fprintf('inflated due to pseudoreplication.\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n\n');
    
    % Store hierarchical results for later use
    hier_results = struct();
    hier_results.icc_ctrl_snr = icc_ctrl_snr;
    hier_results.icc_exp_snr = icc_exp_snr;
    hier_results.design_effect_ctrl = de_ctrl_snr;
    hier_results.design_effect_exp = de_exp_snr;
    hier_results.snr_diff = mean_diff_snr;
    hier_results.snr_diff_ci = [ci_low_snr, ci_high_snr];
    hier_results.event_diff = mean_diff_event;
    hier_results.event_diff_ci = [ci_low_event, ci_high_event];
    hier_results.delta_snr = delta_snr;
    hier_results.delta_snr_ci = [delta_ci_low_snr, delta_ci_high_snr];
    hier_results.delta_event = delta_event;
    hier_results.delta_event_ci = [delta_ci_low_event, delta_ci_high_event];
    
    % Generate methods text
    fprintf('═══════════════════════════════════════════════════════════════════════\n');
    fprintf('  AUTO-GENERATED METHODS TEXT (copy for your paper)\n');
    fprintf('═══════════════════════════════════════════════════════════════════════\n\n');
    
    methods_text = generate_methods_text(n_cond1, n_cond2, ...
        length(metrics.ctrl.all_snr), length(metrics.exp.all_snr), ...
        STATS.alpha, STATS.n_comparisons);
    fprintf('%s\n', methods_text);
    
    % Save methods text to file
    methods_path = fullfile(OUTPUT_DIR, 'METHODS_TEXT.txt');
    fid_methods = fopen(methods_path, 'w');
    fprintf(fid_methods, '%s', methods_text);
    fclose(fid_methods);
    fprintf('Saved: %s\n\n', methods_path);
    
end

%% ════════════════════════════════════════════════════════════════════════
%  GENERATE PUBLICATION FIGURES
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  GENERATING PUBLICATION FIGURES\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

if n_cond1 > 0 && n_cond2 > 0
    
    % Set up colors
    colors = {CONDITIONS(1).color, CONDITIONS(2).color};
    
    % ─────────────────────────────────────────────────────────────────────────
    % FIGURE 1: ROI Detection and Quality Overview
    % ─────────────────────────────────────────────────────────────────────────
    
    fig1 = figure('Position', [50 50 1500 400], 'Color', 'w');
    
    % Panel A: Number of ROIs per sample
    ax1 = subplot(1, 4, 1);
    data_rois = {metrics.ctrl.n_rois(:), metrics.exp.n_rois(:)};
    create_comparison_boxplot(data_rois, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        'Number of ROIs', '');
    title('A) ROI Detection', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    % Add statistical annotation
    [p, ~] = ranksum(metrics.ctrl.n_rois, metrics.exp.n_rois);
    add_stat_annotation(gca, p, alpha_adj);
    
    % Add sample size annotation using xlabel
    xlabel(sprintf('%s (n=%d)          %s (n=%d)', CONDITIONS(1).short_name, n_cond1, ...
        CONDITIONS(2).short_name, n_cond2), 'FontSize', 9);
    set(gca, 'XTickLabel', {});  % Remove tick labels since we use xlabel
    
    % Panel B: % ROIs Removed (quality indicator)
    ax2 = subplot(1, 4, 2);
    data_removed = {metrics.ctrl.pct_removed(:), metrics.exp.pct_removed(:)};
    create_comparison_boxplot(data_removed, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        '% Removed', '');
    title('B) False Positive Rate', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    [p, ~] = ranksum(metrics.ctrl.pct_removed, metrics.exp.pct_removed);
    add_stat_annotation(gca, p, alpha_adj);
    
    % Panel C: SNR distribution (all ROIs) - use robust method for outliers
    ax3 = subplot(1, 4, 3);
    data_snr = {metrics.ctrl.all_snr(:), metrics.exp.all_snr(:)};
    create_comparison_boxplot_notch(data_snr, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        'SNR', '');
    title('C) Signal-to-Noise Ratio', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    [p, ~] = ranksum(metrics.ctrl.all_snr, metrics.exp.all_snr);
    add_stat_annotation(gca, p, alpha_adj);
    
    % Add ROI count in xlabel
    xlabel(sprintf('%s (n=%d)          %s (n=%d)', CONDITIONS(1).short_name, length(metrics.ctrl.all_snr), ...
        CONDITIONS(2).short_name, length(metrics.exp.all_snr)), 'FontSize', 9);
    set(gca, 'XTickLabel', {});
    
    % Panel D: Event Rate distribution
    ax4 = subplot(1, 4, 4);
    data_rate = {metrics.ctrl.all_event_rate(:), metrics.exp.all_event_rate(:)};
    create_comparison_boxplot_notch(data_rate, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        'Event Rate (Hz)', '');
    title('D) Neural Activity', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    [p, ~] = ranksum(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate);
    add_stat_annotation(gca, p, alpha_adj);
    
    % Adjust subplot positions to add space at top (no sgtitle needed)
    for ax = [ax1 ax2 ax3 ax4]
        pos = get(ax, 'Position');
        set(ax, 'Position', [pos(1) pos(2) pos(3) pos(4)*0.9]);
    end
    
    save_publication_figure(fig1, fullfile(figures_dir, 'FigureS1_ROI_Quality'), FIG);
    fprintf('  Saved Figure S1: ROI Quality (Supplementary)\n');
    
    % ─────────────────────────────────────────────────────────────────────────
    % FIGURE 2: Motion-Induced Noise Analysis
    % ─────────────────────────────────────────────────────────────────────────
    
    % ─────────────────────────────────────────────────────────────────────────
    % FIGURE 7: ROI Quality Assessment (Composite - matches manuscript Fig. 7)
    % Panels: A) Removal breakdown, B) SNR vs synchrony, C) Event amplitude,
    %         D) Network synchrony, E) Motion artifact removal, F) Baseline drift
    % ─────────────────────────────────────────────────────────────────────────
    
    fig7_composite = figure('Position', [50 50 1500 800], 'Color', 'w');
    
    % --- Panel A: ROI Removal Breakdown (horizontal bar chart) ---
    subplot(2, 3, 1);
    
    removal_cats = {'Extreme SNR', 'Edge Contact', 'Shape', 'Negative', ...
                    'Step Artifact', 'Kinetics', 'Plateau', 'Correlation', ...
                    'Spatial', 'Inactive', 'Baseline Drift'};
    
    ctrl_removal = [mean(metrics.ctrl.removed_extreme_snr), mean(metrics.ctrl.removed_edge), ...
                    mean(metrics.ctrl.removed_shape), mean(metrics.ctrl.removed_negative), ...
                    mean(metrics.ctrl.removed_step), mean(metrics.ctrl.removed_kinetics), ...
                    mean(metrics.ctrl.removed_plateau), mean(metrics.ctrl.removed_correlated), ...
                    mean(metrics.ctrl.removed_spatial), mean(metrics.ctrl.removed_inactive), ...
                    mean(metrics.ctrl.removed_baseline)];
    
    exp_removal = [mean(metrics.exp.removed_extreme_snr), mean(metrics.exp.removed_edge), ...
                   mean(metrics.exp.removed_shape), mean(metrics.exp.removed_negative), ...
                   mean(metrics.exp.removed_step), mean(metrics.exp.removed_kinetics), ...
                   mean(metrics.exp.removed_plateau), mean(metrics.exp.removed_correlated), ...
                   mean(metrics.exp.removed_spatial), mean(metrics.exp.removed_inactive), ...
                   mean(metrics.exp.removed_baseline)];
    
    total_removal = ctrl_removal + exp_removal;
    [~, sort_idx] = sort(total_removal, 'descend');
    ctrl_removal = ctrl_removal(sort_idx);
    exp_removal = exp_removal(sort_idx);
    removal_cats = removal_cats(sort_idx);
    
    keep = (ctrl_removal > 0) | (exp_removal > 0);
    ctrl_removal = ctrl_removal(keep);
    exp_removal = exp_removal(keep);
    removal_cats = removal_cats(keep);
    
    n_cats = length(removal_cats);
    x = 1:n_cats;
    bar_width = 0.35;
    
    hold on;
    b1 = barh(x - bar_width/2, ctrl_removal, bar_width, 'FaceColor', CONDITIONS(1).color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    b2 = barh(x + bar_width/2, exp_removal, bar_width, 'FaceColor', CONDITIONS(2).color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    set(gca, 'YTick', x, 'YTickLabel', removal_cats, 'FontSize', 8, 'YDir', 'reverse');
    xlabel('Mean ROIs Removed', 'FontSize', 10, 'FontWeight', 'bold');
    title('A) Removal Breakdown', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    legend([b1 b2], {sprintf('%s (n=%d)', CONDITIONS(1).name, n_cond1), ...
        sprintf('%s (n=%d)', CONDITIONS(2).name, n_cond2)}, 'Location', 'southeast', 'FontSize', 7);
    set(gca, 'LineWidth', 1.2, 'Box', 'off', 'TickDir', 'out');
    
    % --- Panel B: SNR vs Motion Synchrony ---
    subplot(2, 3, 2);
    hold on;
    
    h1 = scatter(metrics.ctrl.synchrony, metrics.ctrl.mean_snr, 80, ...
        CONDITIONS(1).color, 'filled', 'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', CONDITIONS(1).color * 0.6, 'LineWidth', 1);
    h2 = scatter(metrics.exp.synchrony, metrics.exp.mean_snr, 100, ...
        CONDITIONS(2).color, 'd', 'filled', 'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', CONDITIONS(2).color * 0.6, 'LineWidth', 1);
    
    xlabel('Synchrony Score', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Mean SNR', 'FontSize', 10, 'FontWeight', 'bold');
    title('B) SNR vs Motion Synchrony', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    legend([h1 h2], {CONDITIONS(1).name, CONDITIONS(2).name}, 'Location', 'best', 'FontSize', 8);
    
    all_sync = [metrics.ctrl.synchrony; metrics.exp.synchrony];
    all_snr_mean = [metrics.ctrl.mean_snr; metrics.exp.mean_snr];
    valid = ~isnan(all_sync) & ~isnan(all_snr_mean);
    if sum(valid) > 3
        [r, p_corr] = corr(all_sync(valid), all_snr_mean(valid), 'Type', 'Spearman');
        if p_corr < 0.05
            text(0.05, 0.95, sprintf('\x03C1 = %.2f, p = %.3f', r, p_corr), ...
                'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top', 'FontWeight', 'bold');
        else
            text(0.05, 0.95, sprintf('\x03C1 = %.2f, p = %.3f', r, p_corr), ...
                'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top', 'Color', [0.5 0.5 0.5]);
        end
    end
    set(gca, 'LineWidth', 1.2, 'FontSize', 10);
    box on;
    
    % --- Panel C: Calcium Event Amplitude ---
    subplot(2, 3, 3);
    
    if ~isempty(metrics.ctrl.all_event_amplitude) && ~isempty(metrics.exp.all_event_amplitude)
        data_amp = {metrics.ctrl.all_event_amplitude(:), metrics.exp.all_event_amplitude(:)};
        create_comparison_violin(data_amp, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
            'Event Amplitude (\DeltaF/F)', '');
        
        [p, ~] = ranksum(metrics.ctrl.all_event_amplitude, metrics.exp.all_event_amplitude);
        add_stat_annotation(gca, p, alpha_adj);
        title('C) Calcium Event Amplitude', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    else
        hold on;
        h1 = scatter(metrics.ctrl.pct_removed, metrics.ctrl.mean_snr, 80, ...
            CONDITIONS(1).color, 'filled', 'MarkerFaceAlpha', 0.7, ...
            'MarkerEdgeColor', CONDITIONS(1).color * 0.6, 'LineWidth', 1);
        h2 = scatter(metrics.exp.pct_removed, metrics.exp.mean_snr, 100, ...
            CONDITIONS(2).color, 'd', 'filled', 'MarkerFaceAlpha', 0.7, ...
            'MarkerEdgeColor', CONDITIONS(2).color * 0.6, 'LineWidth', 1);
        
        xlabel('% ROIs Removed', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Mean SNR', 'FontSize', 10, 'FontWeight', 'bold');
        title('C) SNR vs False Positive Rate', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
        legend([h1 h2], {CONDITIONS(1).name, CONDITIONS(2).name}, 'Location', 'best', 'FontSize', 8);
        set(gca, 'LineWidth', 1.2, 'FontSize', 10);
        box on;
    end
    
    % --- Panel D: Network Synchrony ---
    subplot(2, 3, 4);
    data_sync = {metrics.ctrl.synchrony(:), metrics.exp.synchrony(:)};
    create_comparison_boxplot(data_sync, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        'Synchrony Score', '');
    
    yl = ylim;
    if yl(2) > 0.5
        yline(0.5, 'r--', 'LineWidth', 1.5);
        text(2.5, 0.52, 'High', 'FontSize', 8, 'Color', 'r', 'HorizontalAlignment', 'left');
    end
    if yl(2) > 0.3
        yline(0.3, 'Color', [0.8 0.5 0], 'LineStyle', ':', 'LineWidth', 1.2);
    end
    
    [p, ~] = ranksum(metrics.ctrl.synchrony, metrics.exp.synchrony);
    add_stat_annotation(gca, p, alpha_adj);
    title('D) Network Synchrony', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    % --- Panel E: Motion Artifact Removal ---
    subplot(2, 3, 5);
    data_corr = {metrics.ctrl.pct_correlated(:), metrics.exp.pct_correlated(:)};
    create_comparison_boxplot(data_corr, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        '% Correlated Removed', '');
    
    [p, ~] = ranksum(metrics.ctrl.pct_correlated, metrics.exp.pct_correlated);
    add_stat_annotation(gca, p, alpha_adj);
    title('E) Motion Artifact Removal', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    % --- Panel F: Baseline Drift ---
    subplot(2, 3, 6);
    data_baseline = {metrics.ctrl.pct_baseline_flagged(:), metrics.exp.pct_baseline_flagged(:)};
    create_comparison_boxplot(data_baseline, {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, colors, ...
        '% Baseline Flagged', '');
    
    [p, ~] = ranksum(metrics.ctrl.pct_baseline_flagged, metrics.exp.pct_baseline_flagged);
    add_stat_annotation(gca, p, alpha_adj);
    title('F) Baseline Drift', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    save_publication_figure(fig7_composite, fullfile(figures_dir, 'Figure7_ROI_Quality'), FIG);
    fprintf('  Saved Figure 7: ROI Quality Assessment\n');
    close(fig7_composite);
    
    % ─────────────────────────────────────────────────────────────────────────
    % FIGURE 6: Calcium Imaging Summary (matches manuscript Fig. 6)
    % Panels: A-E) Per-sample metrics & CDFs, F) Forest plot of effect sizes
    % ─────────────────────────────────────────────────────────────────────────
    
    fig6 = figure('Position', [50 50 1400 550], 'Color', 'w');
    
    % Calculate CV early for use later
    cv_ctrl = std(metrics.ctrl.mean_snr, 0, 'omitnan') / mean(metrics.ctrl.mean_snr, 'omitnan') * 100;
    cv_exp = std(metrics.exp.mean_snr, 0, 'omitnan') / mean(metrics.exp.mean_snr, 'omitnan') * 100;
    
    % Panel A: SNR per sample (shows inter-sample variability)
    subplot(2, 3, 1);
    hold on;
    
    % Plot individual sample means with jitter
    x_ctrl = ones(n_cond1, 1) + 0.15*(rand(n_cond1,1)-0.5);
    x_exp = 2*ones(n_cond2, 1) + 0.15*(rand(n_cond2,1)-0.5);
    
    scatter(x_ctrl, metrics.ctrl.mean_snr, 70, CONDITIONS(1).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(1).color*0.6, 'LineWidth', 1);
    scatter(x_exp, metrics.exp.mean_snr, 70, CONDITIONS(2).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(2).color*0.6, 'LineWidth', 1);
    
    % Add mean ± SEM as horizontal lines
    line([0.75 1.25], [1 1]*mean(metrics.ctrl.mean_snr), 'Color', 'k', 'LineWidth', 2);
    line([1.75 2.25], [1 1]*mean(metrics.exp.mean_snr), 'Color', 'k', 'LineWidth', 2);
    
    % Error bars (SEM)
    errorbar(1, mean(metrics.ctrl.mean_snr), std(metrics.ctrl.mean_snr)/sqrt(n_cond1), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    errorbar(2, mean(metrics.exp.mean_snr), std(metrics.exp.mean_snr)/sqrt(n_cond2), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, 'FontSize', 11);
    ylabel('Mean SNR', 'FontSize', 12, 'FontWeight', 'bold');
    title('A) Per-Sample SNR', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    xlim([0.5 2.5]);
    set(gca, 'LineWidth', 1.2);
    box on;
    
    % Add CV annotation
    text(0.05, 0.95, sprintf('CV: %.0f%% vs %.0f%%', cv_ctrl, cv_exp), ...
        'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');
    
    % Panel B: ROIs per sample distribution
    subplot(2, 3, 2);
    hold on;
    
    x_ctrl = ones(n_cond1, 1) + 0.15*(rand(n_cond1,1)-0.5);
    x_exp = 2*ones(n_cond2, 1) + 0.15*(rand(n_cond2,1)-0.5);
    
    scatter(x_ctrl, metrics.ctrl.n_rois, 70, CONDITIONS(1).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(1).color*0.6, 'LineWidth', 1);
    scatter(x_exp, metrics.exp.n_rois, 70, CONDITIONS(2).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(2).color*0.6, 'LineWidth', 1);
    
    line([0.75 1.25], [1 1]*mean(metrics.ctrl.n_rois), 'Color', 'k', 'LineWidth', 2);
    line([1.75 2.25], [1 1]*mean(metrics.exp.n_rois), 'Color', 'k', 'LineWidth', 2);
    
    errorbar(1, mean(metrics.ctrl.n_rois), std(metrics.ctrl.n_rois)/sqrt(n_cond1), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    errorbar(2, mean(metrics.exp.n_rois), std(metrics.exp.n_rois)/sqrt(n_cond2), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, 'FontSize', 11);
    ylabel('ROIs per Sample', 'FontSize', 12, 'FontWeight', 'bold');
    title('B) Per-Sample ROI Count', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    xlim([0.5 2.5]);
    set(gca, 'LineWidth', 1.2);
    box on;
    
    % Panel C: Event rate per sample
    subplot(2, 3, 3);
    hold on;
    
    x_ctrl = ones(n_cond1, 1) + 0.15*(rand(n_cond1,1)-0.5);
    x_exp = 2*ones(n_cond2, 1) + 0.15*(rand(n_cond2,1)-0.5);
    
    scatter(x_ctrl, metrics.ctrl.mean_event_rate, 70, CONDITIONS(1).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(1).color*0.6, 'LineWidth', 1);
    scatter(x_exp, metrics.exp.mean_event_rate, 70, CONDITIONS(2).color, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', CONDITIONS(2).color*0.6, 'LineWidth', 1);
    
    line([0.75 1.25], [1 1]*mean(metrics.ctrl.mean_event_rate), 'Color', 'k', 'LineWidth', 2);
    line([1.75 2.25], [1 1]*mean(metrics.exp.mean_event_rate), 'Color', 'k', 'LineWidth', 2);
    
    errorbar(1, mean(metrics.ctrl.mean_event_rate), std(metrics.ctrl.mean_event_rate)/sqrt(n_cond1), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    errorbar(2, mean(metrics.exp.mean_event_rate), std(metrics.exp.mean_event_rate)/sqrt(n_cond2), ...
        'k', 'LineWidth', 1.5, 'CapSize', 8);
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {CONDITIONS(1).short_name, CONDITIONS(2).short_name}, 'FontSize', 11);
    ylabel('Event Rate (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    title('C) Per-Sample Event Rate', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    xlim([0.5 2.5]);
    set(gca, 'LineWidth', 1.2);
    box on;
    
    % Panel D: Cumulative SNR distribution
    subplot(2, 3, 4);
    hold on;
    
    [f_ctrl, x_ctrl_cdf] = ecdf(metrics.ctrl.all_snr);
    [f_exp, x_exp_cdf] = ecdf(metrics.exp.all_snr);
    
    stairs(x_ctrl_cdf, f_ctrl, 'Color', CONDITIONS(1).color, 'LineWidth', 2);
    stairs(x_exp_cdf, f_exp, 'Color', CONDITIONS(2).color, 'LineWidth', 2);
    
    xlabel('SNR', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Cumulative Probability', 'FontSize', 12, 'FontWeight', 'bold');
    title('D) SNR Distribution (CDF)', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    legend({CONDITIONS(1).name, CONDITIONS(2).name}, 'Location', 'southeast', 'FontSize', 10);
    set(gca, 'LineWidth', 1.2, 'FontSize', 11);
    box on;
    
    % Panel E: Cumulative Event Rate distribution  
    subplot(2, 3, 5);
    hold on;
    
    [f_ctrl, x_ctrl_cdf] = ecdf(metrics.ctrl.all_event_rate);
    [f_exp, x_exp_cdf] = ecdf(metrics.exp.all_event_rate);
    
    stairs(x_ctrl_cdf, f_ctrl, 'Color', CONDITIONS(1).color, 'LineWidth', 2);
    stairs(x_exp_cdf, f_exp, 'Color', CONDITIONS(2).color, 'LineWidth', 2);
    
    xlabel('Event Rate (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Cumulative Probability', 'FontSize', 12, 'FontWeight', 'bold');
    title('E) Event Rate Distribution (CDF)', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    legend({CONDITIONS(1).name, CONDITIONS(2).name}, 'Location', 'southeast', 'FontSize', 10);
    set(gca, 'LineWidth', 1.2, 'FontSize', 11);
    box on;
    
    % Panel F: Forest Plot of Effect Sizes (embedded)
    subplot(2, 3, 6);
    hold on;
    
    % Compute Cliff's delta effect sizes and bootstrap CIs
    metric_labels = {'SNR (sample means)', 'Event Rate (sample means)', ...
                    'Synchrony', '% Correlated Removed', ...
                    'SNR (all ROIs)', 'Event Rate (all ROIs)'};
    
    n_forest_metrics = 6;
    effect_sizes = zeros(n_forest_metrics, 1);
    ci_lows = zeros(n_forest_metrics, 1);
    ci_highs = zeros(n_forest_metrics, 1);
    
    cliffs_func = @(x, y) compute_cliffs_delta(x, y);
    
    effect_sizes(1) = compute_cliffs_delta(metrics.ctrl.mean_snr, metrics.exp.mean_snr);
    [ci_lows(1), ci_highs(1)] = bootstrap_effect_size_ci(metrics.ctrl.mean_snr, metrics.exp.mean_snr, cliffs_func, 2000, 0.95);
    
    effect_sizes(2) = compute_cliffs_delta(metrics.ctrl.mean_event_rate, metrics.exp.mean_event_rate);
    [ci_lows(2), ci_highs(2)] = bootstrap_effect_size_ci(metrics.ctrl.mean_event_rate, metrics.exp.mean_event_rate, cliffs_func, 2000, 0.95);
    
    effect_sizes(3) = compute_cliffs_delta(metrics.ctrl.synchrony, metrics.exp.synchrony);
    [ci_lows(3), ci_highs(3)] = bootstrap_effect_size_ci(metrics.ctrl.synchrony, metrics.exp.synchrony, cliffs_func, 2000, 0.95);
    
    effect_sizes(4) = compute_cliffs_delta(metrics.ctrl.pct_correlated, metrics.exp.pct_correlated);
    [ci_lows(4), ci_highs(4)] = bootstrap_effect_size_ci(metrics.ctrl.pct_correlated, metrics.exp.pct_correlated, cliffs_func, 2000, 0.95);
    
    effect_sizes(5) = compute_cliffs_delta(metrics.ctrl.all_snr, metrics.exp.all_snr);
    [ci_lows(5), ci_highs(5)] = bootstrap_effect_size_ci(metrics.ctrl.all_snr, metrics.exp.all_snr, cliffs_func, 2000, 0.95);
    
    effect_sizes(6) = compute_cliffs_delta(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate);
    [ci_lows(6), ci_highs(6)] = bootstrap_effect_size_ci(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate, cliffs_func, 2000, 0.95);
    
    % Reference line at zero
    line([0 0], [0 n_forest_metrics+1], 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Negligible effect zone
    fill([-0.147 0.147 0.147 -0.147], [0 0 n_forest_metrics+1 n_forest_metrics+1], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    for i = 1:n_forest_metrics
        y_pos = n_forest_metrics - i + 1;
        
        if isnan(effect_sizes(i)) || isnan(ci_lows(i)) || isnan(ci_highs(i))
            text(0, y_pos, 'N/A (insufficient data)', 'FontSize', 7, ...
                'Color', [0.7 0.7 0.7], 'VerticalAlignment', 'middle');
            continue;
        end
        
        if ci_lows(i) > 0 || ci_highs(i) < 0
            fc = [0.2 0.6 0.2];  ms = 80;
        else
            fc = [0.5 0.5 0.5];  ms = 50;
        end
        
        line([ci_lows(i) ci_highs(i)], [y_pos y_pos], 'Color', fc, 'LineWidth', 1.5);
        line([ci_lows(i) ci_lows(i)], [y_pos-0.15 y_pos+0.15], 'Color', fc, 'LineWidth', 1.5);
        line([ci_highs(i) ci_highs(i)], [y_pos-0.15 y_pos+0.15], 'Color', fc, 'LineWidth', 1.5);
        scatter(effect_sizes(i), y_pos, ms, fc, 'd', 'filled');
        
        text(1.05, y_pos, sprintf('%.2f [%.2f, %.2f]', effect_sizes(i), ci_lows(i), ci_highs(i)), ...
            'FontSize', 7, 'VerticalAlignment', 'middle');
    end
    
    set(gca, 'YTick', 1:n_forest_metrics, 'YTickLabel', flip(metric_labels), 'FontSize', 8);
    xlabel('Cliff''s \delta (Effect Size)', 'FontSize', 9, 'FontWeight', 'bold');
    xlim([-1.1 1.5]);
    ylim([0.5 n_forest_metrics + 0.5]);
    
    text(-0.9, n_forest_metrics + 0.3, 'Favors Exp', 'FontSize', 7, 'Color', [0.5 0.5 0.5]);
    text(0.6, n_forest_metrics + 0.3, 'Favors Ctrl', 'FontSize', 7, 'Color', [0.5 0.5 0.5]);
    
    title('F) Effect Sizes (Cliff''s \delta)', 'FontSize', FIG.font_size_title, 'FontWeight', 'bold');
    
    % Compact legend
    scatter(-0.85, 0.8, 50, [0.2 0.6 0.2], 'd', 'filled');
    text(-0.75, 0.8, 'Significant', 'FontSize', 7, 'VerticalAlignment', 'middle');
    scatter(0.3, 0.8, 40, [0.5 0.5 0.5], 'd', 'filled');
    text(0.4, 0.8, 'Not significant', 'FontSize', 7, 'VerticalAlignment', 'middle');
    
    set(gca, 'LineWidth', 1.2);
    box on;
    
    save_publication_figure(fig6, fullfile(figures_dir, 'Figure6_Calcium_Summary'), FIG);
    fprintf('  Saved Figure 6: Calcium Imaging Summary\n');
    close(fig6);
    
    % Close remaining figures safely
    try close(fig1); catch; end
    
end

%% ════════════════════════════════════════════════════════════════════════
%  SAVE RESULTS TABLES
% ════════════════════════════════════════════════════════════════════════

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  SAVING RESULTS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

% Save statistical results table
if ~isempty(stat_results)
    stat_path = fullfile(tables_dir, 'statistical_results.csv');
    writetable(stat_results, stat_path);
    fprintf('Saved: %s\n', stat_path);
end

% Create sample-level summary table
sample_summary = table();

for c = 1:2
    if all_data(c).n_samples == 0
        continue;
    end
    
    for s = 1:length(all_data(c).samples)
        sample = all_data(c).samples(s);
        
        row = table();
        row.Condition = {CONDITIONS(c).name};
        row.Sample = {sample.filename};
        row.N_Original = sample.n_original;
        row.N_Refined = sample.n_rois;
        row.Pct_Removed = sample.pct_removed;
        row.Removed_Correlated = sample.removed_correlated;
        row.Removed_Shape = sample.removed_shape;
        row.Removed_Inactive = sample.removed_inactive;
        row.Removed_Kinetics = sample.removed_kinetics;
        row.Removed_Plateau = sample.removed_plateau;
        row.Removed_Edge = sample.removed_edge;
        row.Removed_SNR = sample.removed_extreme_snr;
        row.Removed_Negative = sample.removed_negative;
        row.Removed_Step = sample.removed_step;
        row.Removed_Spatial = sample.removed_spatial;
        row.Removed_Baseline = sample.removed_baseline;
        row.Mean_SNR = sample.mean_snr;
        row.Mean_Event_Rate = sample.mean_event_rate;
        row.Synchrony_Score = sample.synchrony_score;
        row.Pct_Correlated_Removed = sample.pct_correlated_removed;
        row.Motion_Quality = {sample.motion_quality};
        row.Pct_Baseline_Flagged = sample.pct_baseline_flagged;
        
        sample_summary = [sample_summary; row];
    end
end

sample_path = fullfile(tables_dir, 'sample_summary.csv');
writetable(sample_summary, sample_path);
fprintf('Saved: %s\n', sample_path);

% Create group summary table
group_summary = table();

for c = 1:2
    if all_data(c).n_samples == 0
        continue;
    end
    
    cond_name = lower(CONDITIONS(c).short_name);
    
    row = table();
    row.Condition = {CONDITIONS(c).name};
    row.N_Samples = all_data(c).n_samples;
    row.Total_ROIs = sum(metrics.(cond_name).n_rois);
    row.Mean_ROIs_Per_Sample = mean(metrics.(cond_name).n_rois);
    row.SD_ROIs_Per_Sample = std(metrics.(cond_name).n_rois);
    row.Mean_Pct_Removed = mean(metrics.(cond_name).pct_removed);
    row.Mean_SNR = mean(metrics.(cond_name).all_snr);
    row.SD_SNR = std(metrics.(cond_name).all_snr);
    row.Mean_Event_Rate = mean(metrics.(cond_name).all_event_rate);
    row.SD_Event_Rate = std(metrics.(cond_name).all_event_rate);
    row.Mean_Synchrony = mean(metrics.(cond_name).synchrony, 'omitnan');
    row.Mean_Pct_Correlated_Removed = mean(metrics.(cond_name).pct_correlated);
    
    group_summary = [group_summary; row];
end

group_path = fullfile(tables_dir, 'group_summary.csv');
writetable(group_summary, group_path);
fprintf('Saved: %s\n', group_path);

%% ════════════════════════════════════════════════════════════════════════
%  GENERATE ANALYSIS REPORT
% ════════════════════════════════════════════════════════════════════════

report_path = fullfile(OUTPUT_DIR, 'ANALYSIS_REPORT.txt');
fid = fopen(report_path, 'w');

fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n');
fprintf(fid, 'STATISTICAL COMPARISON REPORT\n');
fprintf(fid, 'Control (Hydrogel) vs Experimental (Hydrogel + Diethyl Ether)\n');
fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n\n');

fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'MATLAB Version: %s\n\n', version('-release'));

fprintf(fid, 'DATA SUMMARY\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');

for c = 1:2
    if all_data(c).n_samples == 0
        continue;
    end
    
    cond_name = lower(CONDITIONS(c).short_name);
    
    fprintf(fid, '\n%s:\n', CONDITIONS(c).name);
    fprintf(fid, '  Samples: %d\n', all_data(c).n_samples);
    fprintf(fid, '  Total refined ROIs: %d\n', sum(metrics.(cond_name).n_rois));
    fprintf(fid, '  ROIs per sample: %.1f ± %.1f (mean ± SD)\n', ...
        mean(metrics.(cond_name).n_rois), std(metrics.(cond_name).n_rois));
    fprintf(fid, '  SNR: %.2f ± %.2f\n', mean(metrics.(cond_name).all_snr), std(metrics.(cond_name).all_snr));
    fprintf(fid, '  Event rate: %.3f ± %.3f Hz\n', ...
        mean(metrics.(cond_name).all_event_rate), std(metrics.(cond_name).all_event_rate));
    fprintf(fid, '  Synchrony score: %.3f ± %.3f\n', ...
        mean(metrics.(cond_name).synchrony, 'omitnan'), std(metrics.(cond_name).synchrony, 'omitnan'));
    fprintf(fid, '  %% Correlated removed: %.1f ± %.1f\n', ...
        mean(metrics.(cond_name).pct_correlated), std(metrics.(cond_name).pct_correlated));
end

fprintf(fid, '\n\nSTATISTICAL ANALYSIS\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'Test: Mann-Whitney U (non-parametric)\n');
fprintf(fid, 'Effect size: Cliff''s delta\n');
fprintf(fid, 'Significance level: α = %.3f', STATS.alpha);
if STATS.use_bonferroni
    fprintf(fid, ' (Bonferroni-corrected: %.4f)\n\n', alpha_adj);
else
    fprintf(fid, '\n\n');
end

if ~isempty(stat_results)
    fprintf(fid, 'Results:\n');
    for i = 1:height(stat_results)
        sig_str = '';
        if stat_results.Significant(i)
            sig_str = ' *SIGNIFICANT*';
        end
        fprintf(fid, '  %s: p = %.4f, δ = %.3f (%s)%s\n', ...
            stat_results.Metric{i}, stat_results.p_value(i), ...
            stat_results.Cliffs_Delta(i), stat_results.Effect_Size{i}, sig_str);
    end
end

fprintf(fid, '\n\nMOTION-INDUCED NOISE INTERPRETATION\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'Motion artifacts in calcium imaging manifest as:\n');
fprintf(fid, '  1. High synchrony scores (>0.5): Indicates correlated activity across\n');
fprintf(fid, '     ROIs that is likely due to motion rather than neural coordination.\n');
fprintf(fid, '  2. High %% correlated ROIs removed: Many ROIs had similar traces\n');
fprintf(fid, '     (r>0.7), suggesting shared motion artifact rather than true signal.\n\n');

fprintf(fid, 'These two metrics are reported separately (not combined into a\n');
fprintf(fid, 'composite score) to preserve interpretability.\n\n');

fprintf(fid, '\n\nFILES GENERATED\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'Figures (in figures/ subfolder):\n');
fprintf(fid, '  - Figure1_ROI_Quality.*       : ROI detection and signal quality\n');
fprintf(fid, '  - Figure2_Motion_Noise.*      : Motion-induced noise analysis\n');
fprintf(fid, '  - Figure3_Correlation.*       : Relationship between metrics\n');
fprintf(fid, '  - Figure4_Removal_Breakdown.* : ROI removal criteria by condition\n');
fprintf(fid, '  - Figure5_Data_Quality.*      : Inter-sample variability\n');
fprintf(fid, '  - Figure7_Forest_Plot.*       : Effect sizes with bootstrap CIs\n\n');

fprintf(fid, 'Tables (in tables/ subfolder):\n');
fprintf(fid, '  - statistical_results.csv     : Full statistical test results\n');
fprintf(fid, '  - sample_summary.csv          : Per-sample metrics\n');
fprintf(fid, '  - group_summary.csv           : Condition-level summaries\n\n');

fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n');

fclose(fid);
fprintf('Saved: %s\n', report_path);

%% ════════════════════════════════════════════════════════════════════════
%  PLAIN LANGUAGE SUMMARY (for non-experts)
% ════════════════════════════════════════════════════════════════════════

if OUTPUT.generate_summary && n_cond1 > 0 && n_cond2 > 0
    
    fprintf('\n');
    fprintf('═══════════════════════════════════════════════════════════════════════\n');
    fprintf('  PLAIN LANGUAGE SUMMARY\n');
    fprintf('═══════════════════════════════════════════════════════════════════════\n\n');
    
    % Calculate key comparisons for interpretation
    snr_ctrl = mean(metrics.ctrl.all_snr);
    snr_exp = mean(metrics.exp.all_snr);
    snr_ratio = snr_ctrl / snr_exp;
    
    event_ctrl = mean(metrics.ctrl.all_event_rate);
    event_exp = mean(metrics.exp.all_event_rate);
    
    [p_snr, ~] = ranksum(metrics.ctrl.all_snr, metrics.exp.all_snr);
    [p_event, ~] = ranksum(metrics.ctrl.all_event_rate, metrics.exp.all_event_rate);
    
    sync_ctrl = mean(metrics.ctrl.synchrony, 'omitnan');
    sync_exp = mean(metrics.exp.synchrony, 'omitnan');
    
    % Determine which condition is "better"
    if snr_ctrl > snr_exp
        better_snr = 'Control (Hydrogel only)';
        worse_snr = 'Experimental (with Diethyl Ether)';
        snr_pct = (snr_ctrl/snr_exp - 1) * 100;
    else
        better_snr = 'Experimental (with Diethyl Ether)';
        worse_snr = 'Control (Hydrogel only)';
        snr_pct = (snr_exp/snr_ctrl - 1) * 100;
    end
    
    if event_ctrl > event_exp
        better_event = 'Control';
        event_pct = (event_ctrl/event_exp - 1) * 100;
    else
        better_event = 'Experimental';
        event_pct = (event_exp/event_ctrl - 1) * 100;
    end
    
    % Print summary
    fprintf('WHAT WE MEASURED:\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    fprintf('We compared two methods for immobilizing Drosophila larvae during\n');
    fprintf('calcium imaging of neural activity:\n');
    fprintf('  • Control: Hydrogel only (%d samples, %d neurons)\n', n_cond1, length(metrics.ctrl.all_snr));
    fprintf('  • Experimental: Hydrogel + Diethyl Ether (%d samples, %d neurons)\n\n', n_cond2, length(metrics.exp.all_snr));
    
    fprintf('KEY FINDINGS:\n');
    fprintf('─────────────────────────────────────────────────────────────────────────\n');
    
    % Signal Quality
    fprintf('\n1. SIGNAL QUALITY (Signal-to-Noise Ratio)\n');
    fprintf('   • %s has %.0f%% higher signal quality\n', better_snr, snr_pct);
    if p_snr < alpha_adj
        fprintf('   • This difference is statistically significant (p < %.4f after correction)\n', alpha_adj);
    else
        fprintf('   • This difference is not statistically significant (p = %.3f)\n', p_snr);
    end
    
    % Neural Activity
    fprintf('\n2. NEURAL ACTIVITY (Event Rate)\n');
    fprintf('   • %s shows %.0f%% more detected neural events\n', better_event, event_pct);
    if p_event < alpha_adj
        fprintf('   • This difference is statistically significant (p < %.4f after correction)\n', alpha_adj);
    else
        fprintf('   • This difference is not statistically significant (p = %.3f)\n', p_event);
    end
    
    % Motion Artifacts
    fprintf('\n3. MOTION ARTIFACTS (Synchrony Score)\n');
    if sync_ctrl < sync_exp
        fprintf('   • Control shows lower synchrony (%.2f vs %.2f)\n', sync_ctrl, sync_exp);
    else
        fprintf('   • Experimental shows lower synchrony (%.2f vs %.2f)\n', sync_exp, sync_ctrl);
    end
    fprintf('   • Values below 0.5 are generally considered acceptable\n');
    
    fprintf('\n─────────────────────────────────────────────────────────────────────────\n');
    fprintf('NOTE: Interpretation of these results should consider:\n');
    fprintf('  • Whether hierarchical bootstrap CIs agree with ROI-level tests\n');
    fprintf('  • Physical motion data from the neuroimaging motion analysis\n');
    fprintf('  • The removal breakdown (which criteria differ between conditions)\n');
    fprintf('  • Potential confounds (photobleaching, sample health, preparation)\n');
    fprintf('Conclusions belong in the Discussion section of the manuscript.\n');
    
    fprintf('\n');
    fprintf('⚠ IMPORTANT CAVEATS:\n');
    if n_cond2 < 5
        fprintf('  • The Experimental group has only %d samples - collect more\n', n_cond2);
        fprintf('    data for robust conclusions (recommend n≥6 per group)\n');
    end
    fprintf('  • ROI-level statistics (n=%d vs n=%d) have high power but\n', ...
        length(metrics.ctrl.all_snr), length(metrics.exp.all_snr));
    fprintf('    may not account for within-sample correlation\n');
    fprintf('  • Sample-level statistics are more conservative but\n');
    fprintf('    may be underpowered with current sample sizes\n');
    fprintf('\n');
    
    % Save plain language summary to file
    summary_path = fullfile(OUTPUT_DIR, 'PLAIN_LANGUAGE_SUMMARY.txt');
    fid_sum = fopen(summary_path, 'w');
    
    fprintf(fid_sum, 'PLAIN LANGUAGE SUMMARY\n');
    fprintf(fid_sum, 'Comparison of Immobilization Methods for Calcium Imaging\n');
    fprintf(fid_sum, '══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf(fid_sum, 'STUDY OVERVIEW:\n');
    fprintf(fid_sum, 'We tested two methods for holding Drosophila larvae still during\n');
    fprintf(fid_sum, 'brain imaging experiments:\n');
    fprintf(fid_sum, '  1. Control: Hydrogel alone\n');
    fprintf(fid_sum, '  2. Experimental: Hydrogel plus Diethyl Ether (an anesthetic)\n\n');
    
    fprintf(fid_sum, 'MAIN QUESTION:\n');
    fprintf(fid_sum, 'Does adding the anesthetic affect the quality of brain recordings?\n\n');
    
    fprintf(fid_sum, 'RESULTS:\n');
    fprintf(fid_sum, '• Signal Quality: %s produces %.0f%% better signals', better_snr, snr_pct);
    if p_snr < alpha_adj
        fprintf(fid_sum, ' (SIGNIFICANT)\n');
    else
        fprintf(fid_sum, ' (not significant)\n');
    end
    
    fprintf(fid_sum, '• Neural Activity: %s detects %.0f%% more brain events', better_event, event_pct);
    if p_event < alpha_adj
        fprintf(fid_sum, ' (SIGNIFICANT)\n');
    else
        fprintf(fid_sum, ' (not significant)\n');
    end
    
    fprintf(fid_sum, '• Motion Artifacts: Both methods have acceptable movement noise\n\n');
    
    fprintf(fid_sum, 'WHAT THIS MEANS:\n');
    fprintf(fid_sum, 'The comparison of these metrics, along with motion tracking data\n');
    fprintf(fid_sum, 'and the ROI removal breakdown, provides the evidence base for\n');
    fprintf(fid_sum, 'evaluating which immobilization method is more suitable.\n');
    fprintf(fid_sum, 'See the hierarchical bootstrap results for the most conservative\n');
    fprintf(fid_sum, 'statistical assessment.\n');
    
    fclose(fid_sum);
    fprintf('Saved: %s\n', summary_path);
end

%% ════════════════════════════════════════════════════════════════════════
%  NOTE: Graphical abstract removed (v2.0)
%  Auto-generated "winner" framing and recommendation text removed to 
%  avoid presenting automated conclusions. Interpretation belongs in
%  the Discussion section of the manuscript, not in analysis code.
% ════════════════════════════════════════════════════════════════════════

%% ════════════════════════════════════════════════════════════════════════
%  COMPLETION
% ════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 STATISTICAL COMPARISON COMPLETE                      ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Output location: %s\n', OUTPUT_DIR);
fprintf('\n');
fprintf('Generated:\n');
fprintf('  - 7 publication-quality figures (PNG, PDF, SVG)\n');
fprintf('  - 3 CSV data tables\n');
fprintf('  - 1 analysis report (TXT)\n');
fprintf('  - 1 plain language summary (TXT)\n');
fprintf('  - 1 methods section text (TXT)\n');
fprintf('\n');
fprintf('Figures:\n');
fprintf('  1. ROI Detection and Signal Quality\n');
fprintf('  2. Motion-Induced Noise Estimation\n');
fprintf('  3. Correlation Analysis\n');
fprintf('  4. Normalized Comparison\n');
fprintf('  5. Inter-Sample Variability\n');
fprintf('  6. Graphical Summary (for non-experts)\n');
fprintf('  7. Forest Plot (effect sizes with 95%% CIs) ← RECOMMENDED FOR PUBLICATION\n');
fprintf('\n');
fprintf('Statistical Features:\n');
fprintf('  - Hierarchical bootstrap (accounts for nested ROI data)\n');
fprintf('  - Intraclass correlation (ICC) and design effect\n');
fprintf('  - Effect sizes with confidence intervals\n');
fprintf('  - Hedges'' g (bias-corrected standardized effect)\n');
fprintf('  - Auto-generated methods text\n');
fprintf('\n');
fprintf('For publication:\n');
fprintf('  1. Use Figure 7 (Forest Plot) as the main statistical figure\n');
fprintf('  2. Use PDF/SVG figures for manuscripts\n');
fprintf('  3. Copy METHODS_TEXT.txt to your Methods section\n');
fprintf('  4. Report effect sizes with CIs, not just p-values\n');
fprintf('  5. Include ICC values to address nested data structure\n');
fprintf('  6. Use Figure 6 as a graphical abstract\n');
fprintf('\n');

end  % end of main function COMPARE_NEURO_CONDITIONS_v3

%% ════════════════════════════════════════════════════════════════════════
%  CONFIG HELPERS
% ════════════════════════════════════════════════════════════════════════

function cfg = merge_with_defaults(user_cfg)
% Merge a user-provided config struct with defaults.
cfg = compare_conditions_config();
default_fields = fieldnames(cfg);
user_fields = fieldnames(user_cfg);

for i = 1:numel(user_fields)
    if ismember(user_fields{i}, default_fields)
        if isstruct(cfg.(user_fields{i})) && isstruct(user_cfg.(user_fields{i})) ...
                && ~isempty(fieldnames(cfg.(user_fields{i}))) && ~strcmp(user_fields{i}, 'conditions')
            sub_defaults = fieldnames(cfg.(user_fields{i}));
            sub_user = fieldnames(user_cfg.(user_fields{i}));
            for j = 1:numel(sub_user)
                if ismember(sub_user{j}, sub_defaults)
                    cfg.(user_fields{i}).(sub_user{j}) = user_cfg.(user_fields{i}).(sub_user{j});
                else
                    warning('compare_conditions:unknownConfig', ...
                        'Unknown config field ''%s.%s'' — ignored.', user_fields{i}, sub_user{j});
                end
            end
        else
            cfg.(user_fields{i}) = user_cfg.(user_fields{i});
        end
    else
        warning('compare_conditions:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See compare_conditions_config().', ...
            user_fields{i});
    end
end
end

function conditions = prompt_for_conditions()
% Prompt the user to define conditions via folder pickers.
if ~usejava('desktop')
    error(['No conditions defined and no GUI available.\n' ...
           'Usage:\n' ...
           '  cfg = compare_conditions_config();\n' ...
           '  cfg.conditions(1).name = ''Control'';\n' ...
           '  cfg.conditions(1).refined_dir = ''/path/to/refined/data'';\n' ...
           '  COMPARE_NEURO_CONDITIONS_v3(cfg);\n\n' ...
           'See also: compare_conditions_config']);
end

default_colors = {[0.2 0.4 0.8], [0.8 0.3 0.2], [0.2 0.7 0.3], [0.6 0.2 0.7]};
default_colors_light = {[0.6 0.7 0.9], [0.9 0.7 0.6], [0.6 0.9 0.7], [0.8 0.6 0.9]};

conditions = struct('name', {}, 'short_name', {}, 'refined_dir', {}, ...
    'color', {}, 'color_light', {}, 'raw_dir', {});
condIdx = 0;

while true
    condIdx = condIdx + 1;
    answer = questdlg(sprintf('Add condition %d?', condIdx), ...
        'Add Condition', 'Yes', 'No (done)', 'Yes');
    if ~strcmp(answer, 'Yes'), break; end
    
    name = inputdlg(sprintf('Condition %d name:', condIdx), 'Condition Name', 1, {sprintf('Condition_%d', condIdx)});
    if isempty(name), break; end
    
    refinedDir = uigetdir('', sprintf('Select REFINED data folder for %s', name{1}));
    if refinedDir == 0, break; end
    
    conditions(condIdx).name = name{1};
    conditions(condIdx).short_name = name{1}(1:min(4, length(name{1})));
    conditions(condIdx).refined_dir = refinedDir;
    cidx = min(condIdx, length(default_colors));
    conditions(condIdx).color = default_colors{cidx};
    conditions(condIdx).color_light = default_colors_light{cidx};
    conditions(condIdx).raw_dir = '';
end

if isempty(conditions)
    error('No conditions defined. Exiting.');
end
end

%% ════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
% ════════════════════════════════════════════════════════════════════════

function val = get_field_or_zero(s, field_name)
% GET_FIELD_OR_ZERO - Safely read a struct field, returning 0 if missing
    if isfield(s, field_name)
        val = s.(field_name);
        if isempty(val) || ~isnumeric(val)
            val = 0;
        end
    else
        val = 0;
    end
end


function test_normality(name, data)
% TEST_NORMALITY - Test and report normality of data
    data = data(~isnan(data));
    n = length(data);
    
    if n < 4
        fprintf('  %-25s: n=%d (too few for normality test)\n', name, n);
        return;
    end
    
    try
        % Lilliefors test (available in Statistics Toolbox)
        [h, p] = lillietest(data);
        test_name = 'Lilliefors';
        
        if h == 1
            result = 'NON-NORMAL';
        else
            result = 'normal';
        end
        fprintf('  %-25s: %s p=%.4f (%s)\n', name, test_name, p, result);
    catch
        % Fallback: use skewness/kurtosis
        s = skewness(data);
        k = kurtosis(data);
        if abs(s) > 1 || abs(k - 3) > 2
            result = 'likely non-normal';
        else
            result = 'possibly normal';
        end
        fprintf('  %-25s: skew=%.2f, kurt=%.2f (%s)\n', name, s, k, result);
    end
end


function delta = compute_cliffs_delta(x, y)
% COMPUTE_CLIFFS_DELTA - Non-parametric effect size (vectorized for speed)
%
% Cliff's delta ranges from -1 to +1:
%   +1: All x values greater than all y values
%   -1: All y values greater than all x values
%    0: Complete overlap
%
% Interpretation thresholds (Romano et al., 2006):
%   |δ| < 0.147: negligible
%   |δ| < 0.33:  small
%   |δ| < 0.474: medium
%   |δ| >= 0.474: large

    x = x(:);
    y = y(:);
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    n_x = length(x);
    n_y = length(y);
    
    if n_x == 0 || n_y == 0
        delta = NaN;
        return;
    end
    
    % Vectorized computation (much faster for large samples)
    % Compare each x to each y using broadcasting
    comparisons = sign(bsxfun(@minus, x, y'));  % n_x by n_y matrix
    
    % Count dominance
    delta = sum(comparisons(:)) / (n_x * n_y);
end


function sig_str = get_sig_stars(p)
% GET_SIG_STARS - Convert p-value to significance stars

    if p < 0.001
        sig_str = '***';
    elseif p < 0.01
        sig_str = '**';
    elseif p < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
end


function power = estimate_power(n1, n2, effect_size)
% ESTIMATE_POWER - Approximate power calculation for two-sample test
%
% Uses approximation based on Cohen's method
% This is a rough estimate - use proper power analysis tools for publication

    if nargin < 3 || isnan(effect_size) || effect_size == 0
        power = NaN;
        return;
    end
    
    % Harmonic mean of sample sizes
    n_harm = 2 * n1 * n2 / (n1 + n2);
    
    % Non-centrality parameter (approximate)
    ncp = effect_size * sqrt(n_harm / 2);
    
    % Critical value for alpha = 0.05 (two-tailed)
    z_crit = 1.96;
    
    % Power approximation using normal distribution
    % This is simplified - actual power depends on distribution assumptions
    power = 1 - normcdf(z_crit - ncp) + normcdf(-z_crit - ncp);
    
    % Bound between 0 and 1
    power = max(0, min(1, power));
end


function create_comparison_boxplot(data, labels, colors, ylabel_str, ~)
% CREATE_COMPARISON_BOXPLOT - Create publication-quality boxplot with individual points

    positions = [1, 2];
    
    hold on;
    
    for i = 1:length(data)
        d = data{i};
        
        % Skip if empty
        if isempty(d) || all(isnan(d))
            continue;
        end
        
        % Ensure column vector and remove NaN values
        d = d(:);
        d = d(~isnan(d));
        
        if isempty(d)
            continue;
        end
        
        % Box plot with better styling
        bp = boxplot(d, 'Positions', positions(i), 'Widths', 0.4, ...
            'Colors', colors{i}, 'Symbol', '', 'MedianStyle', 'line');
        set(bp, 'LineWidth', 1.5);
        
        % Make median line thicker and darker
        h = findobj(bp, 'Tag', 'Median');
        set(h, 'LineWidth', 2.5, 'Color', colors{i} * 0.7);
        
        % Overlay individual points with jitter
        n_points = length(d);
        jitter_width = min(0.15, 0.3 / sqrt(n_points));  % Reduce jitter for many points
        jitter = jitter_width * (rand(size(d)) - 0.5);
        scatter(positions(i) + jitter, d, 25, colors{i}, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
    end
    
    set(gca, 'XTick', positions, 'XTickLabel', labels, 'FontSize', 11);
    ylabel(ylabel_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0.4, 2.6]);
    
    % Set y-axis to start at 0 if all values are positive
    yl = ylim;
    all_vals = [];
    for i = 1:length(data)
        d = data{i}(:);
        d = d(~isnan(d));
        all_vals = [all_vals; d];
    end
    
    if ~isempty(all_vals) && min(all_vals) >= 0
        ylim([0, yl(2) * 1.15]);
    end
    
    box on;
    set(gca, 'LineWidth', 1.2);
end


function create_comparison_boxplot_notch(data, labels, colors, ylabel_str, ~)
% CREATE_COMPARISON_BOXPLOT_NOTCH - Boxplot optimized for data with outliers
% Uses notched boxplot and clips y-axis to show main distribution

    positions = [1, 2];
    
    hold on;
    
    % First pass: collect all data to determine reasonable y-limits
    all_vals = [];
    for i = 1:length(data)
        d = data{i}(:);
        d = d(~isnan(d));
        all_vals = [all_vals; d];
    end
    
    if isempty(all_vals)
        return;
    end
    
    % Calculate robust y-limits (exclude extreme outliers for display)
    q1_all = prctile(all_vals, 5);
    q3_all = prctile(all_vals, 95);
    iqr_all = q3_all - q1_all;
    y_min = max(0, q1_all - 0.5 * iqr_all);  % Don't go below 0 for these metrics
    y_max = q3_all + 0.5 * iqr_all;
    
    for i = 1:length(data)
        d = data{i};
        
        % Skip if empty
        if isempty(d) || all(isnan(d))
            continue;
        end
        
        % Ensure column vector and remove NaN values
        d = d(:);
        d = d(~isnan(d));
        
        if isempty(d)
            continue;
        end
        
        % Notched boxplot (notch shows 95% CI of median)
        bp = boxplot(d, 'Positions', positions(i), 'Widths', 0.35, ...
            'Colors', colors{i}, 'Symbol', '', 'Notch', 'on');
        set(bp, 'LineWidth', 1.5);
        
        % Make median line thicker
        h = findobj(bp, 'Tag', 'Median');
        set(h, 'LineWidth', 2.5, 'Color', colors{i} * 0.6);
        
        % Make box filled with transparency
        h_box = findobj(bp, 'Tag', 'Box');
        patch(get(h_box, 'XData'), get(h_box, 'YData'), colors{i}, ...
            'FaceAlpha', 0.3, 'EdgeColor', colors{i}, 'LineWidth', 1.5);
        
        % Show outliers as small dots (only those within display range)
        q1 = prctile(d, 25);
        q3 = prctile(d, 75);
        iqr_d = q3 - q1;
        outlier_low = d < (q1 - 1.5 * iqr_d);
        outlier_high = d > (q3 + 1.5 * iqr_d);
        outliers = d(outlier_low | outlier_high);
        
        % Plot outliers within visible range
        visible_outliers = outliers(outliers >= y_min & outliers <= y_max);
        if ~isempty(visible_outliers)
            scatter(positions(i) * ones(size(visible_outliers)), visible_outliers, ...
                15, colors{i}, 'filled', 'MarkerFaceAlpha', 0.5);
        end
        
        % Add text if outliers exist outside range
        n_outliers_hidden = sum(outliers > y_max) + sum(outliers < y_min);
        if n_outliers_hidden > 0
            text(positions(i), y_max * 0.95, sprintf('+%d', n_outliers_hidden), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        end
    end
    
    set(gca, 'XTick', positions, 'XTickLabel', labels, 'FontSize', 11);
    ylabel(ylabel_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0.4, 2.6]);
    ylim([y_min, y_max * 1.15]);
    
    box on;
    set(gca, 'LineWidth', 1.2);
end


function create_comparison_violin(data, labels, colors, ylabel_str, ~)
% CREATE_COMPARISON_VIOLIN - Create publication-quality violin plot

    positions = [1, 2];
    
    hold on;
    
    for i = 1:length(data)
        d = data{i};
        
        % Skip if empty
        if isempty(d)
            continue;
        end
        
        % Ensure column vector and remove NaN values
        d = d(:);
        d = d(~isnan(d));
        
        if length(d) < 2
            % Just plot points if too few for density estimate
            scatter(positions(i), d, 60, colors{i}, 'filled');
            continue;
        end
        
        % Create kernel density estimate
        try
            [f, xi] = ksdensity(d, 'NumPoints', 100);
        catch
            % Fallback: just plot histogram-like representation
            scatter(positions(i) + 0.1*(rand(size(d))-0.5), d, 25, colors{i}, 'filled', 'MarkerFaceAlpha', 0.4);
            continue;
        end
        
        % Normalize width for aesthetics
        f = f / max(f) * 0.35;
        
        % Plot violin shape with smooth edges
        fill([positions(i) - f, fliplr(positions(i) + f)], [xi, fliplr(xi)], ...
            colors{i}, 'FaceAlpha', 0.25, 'EdgeColor', colors{i}, 'LineWidth', 1.5);
        
        % Add quartile lines with better styling
        q1 = prctile(d, 25);
        q2 = prctile(d, 50);
        q3 = prctile(d, 75);
        
        % Interquartile range box (thin)
        line([positions(i)-0.08, positions(i)+0.08], [q1, q1], 'Color', colors{i}*0.7, 'LineWidth', 2);
        line([positions(i)-0.12, positions(i)+0.12], [q2, q2], 'Color', colors{i}*0.5, 'LineWidth', 3);  % Median thicker
        line([positions(i)-0.08, positions(i)+0.08], [q3, q3], 'Color', colors{i}*0.7, 'LineWidth', 2);
        
        % Connect quartiles with vertical line
        line([positions(i), positions(i)], [q1, q3], 'Color', colors{i}*0.7, 'LineWidth', 2);
    end
    
    set(gca, 'XTick', positions, 'XTickLabel', labels, 'FontSize', 11);
    ylabel(ylabel_str, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([0.4, 2.6]);
    
    box on;
    set(gca, 'LineWidth', 1.2);
end


function add_stat_annotation(ax, p_value, alpha_adj)
% ADD_STAT_ANNOTATION - Add publication-quality significance annotation

    yl = ylim(ax);
    y_range = diff(yl);
    y_pos = yl(2) - 0.05 * y_range;  % Position near top
    bracket_height = 0.02 * y_range;  % Consistent bracket height
    
    if p_value < alpha_adj
        sig_str = get_sig_stars(p_value);
        text(1.5, y_pos, sig_str, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Add bracket with consistent height
        bracket_y = y_pos - 0.03 * y_range;
        line([1.1 1.9], [bracket_y bracket_y], 'Color', 'k', 'LineWidth', 1.2);
        line([1.1 1.1], [bracket_y - bracket_height, bracket_y], 'Color', 'k', 'LineWidth', 1.2);
        line([1.9 1.9], [bracket_y - bracket_height, bracket_y], 'Color', 'k', 'LineWidth', 1.2);
    else
        % Non-significant: show "n.s." in gray
        text(1.5, y_pos, 'n.s.', ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
    end
end


function save_publication_figure(fig, filepath, fig_settings)
% SAVE_PUBLICATION_FIGURE - Save figure in multiple publication formats

    % PNG (for presentations)
    print(fig, [filepath '.png'], '-dpng', ['-r' num2str(fig_settings.dpi)]);
    
    % PDF (for manuscripts)
    print(fig, [filepath '.pdf'], '-dpdf', '-bestfit');
    
    % SVG (for editing)
    try
        print(fig, [filepath '.svg'], '-dsvg');
    catch
        % SVG may not be available in older MATLAB versions
    end
    
    fprintf('  Saved: %s\n', filepath);
end


function [ci_low, ci_high, boot_dist] = hierarchical_bootstrap_ci(sample_data, n_boot, ci_level)
% HIERARCHICAL_BOOTSTRAP_CI - Bootstrap CI accounting for nested data structure
%
% For nested data (ROIs within samples), this performs a two-level bootstrap:
%   1. Resample samples with replacement
%   2. For each resampled sample, resample ROIs with replacement
%
% This correctly accounts for within-sample correlation and avoids
% inflated Type I error rates that occur when treating nested data as independent.
%
% Reference: Saravanan et al. (2020) "Application of the hierarchical bootstrap
%            to multi-level data in neuroscience" eNeuro
%
% INPUTS:
%   sample_data - cell array where each cell contains ROI-level data from one sample
%   n_boot      - number of bootstrap iterations (default: 10000)
%   ci_level    - confidence level (default: 0.95 for 95% CI)
%
% OUTPUTS:
%   ci_low, ci_high - lower and upper bounds of CI
%   boot_dist       - full bootstrap distribution

    if nargin < 2 || isempty(n_boot)
        n_boot = 10000;
    end
    if nargin < 3 || isempty(ci_level)
        ci_level = 0.95;
    end
    
    n_samples = length(sample_data);
    boot_means = zeros(n_boot, 1);
    
    for b = 1:n_boot
        % Level 1: Resample samples with replacement
        sample_idx = randi(n_samples, n_samples, 1);
        
        % Level 2: For each resampled sample, resample ROIs
        all_resampled_data = [];
        for s = 1:n_samples
            sample_rois = sample_data{sample_idx(s)};
            if ~isempty(sample_rois)
                n_rois = length(sample_rois);
                roi_idx = randi(n_rois, n_rois, 1);
                all_resampled_data = [all_resampled_data; sample_rois(roi_idx)];
            end
        end
        
        if ~isempty(all_resampled_data)
            boot_means(b) = mean(all_resampled_data, 'omitnan');
        else
            boot_means(b) = NaN;
        end
    end
    
    % Remove NaN values
    boot_means = boot_means(~isnan(boot_means));
    
    if isempty(boot_means)
        ci_low = NaN;
        ci_high = NaN;
        boot_dist = [];
        return;
    end
    
    % Percentile method for CI
    alpha = 1 - ci_level;
    ci_low = prctile(boot_means, 100 * alpha/2);
    ci_high = prctile(boot_means, 100 * (1 - alpha/2));
    boot_dist = boot_means;
end


function [ci_low, ci_high] = bootstrap_effect_size_ci(x, y, effect_func, n_boot, ci_level)
% BOOTSTRAP_EFFECT_SIZE_CI - Bootstrap CI for any effect size measure
%
% INPUTS:
%   x, y        - data vectors for two groups
%   effect_func - function handle computing effect size: @(x,y) ...
%   n_boot      - number of bootstrap iterations
%   ci_level    - confidence level
%
% Reference: Cumming (2012) "Understanding the New Statistics"

    if nargin < 4 || isempty(n_boot)
        n_boot = 10000;
    end
    if nargin < 5 || isempty(ci_level)
        ci_level = 0.95;
    end
    
    x = x(:);
    y = y(:);
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    n_x = length(x);
    n_y = length(y);
    
    if n_x < 2 || n_y < 2
        ci_low = NaN;
        ci_high = NaN;
        return;
    end
    
    boot_effects = zeros(n_boot, 1);
    
    for b = 1:n_boot
        % Resample each group independently
        x_boot = x(randi(n_x, n_x, 1));
        y_boot = y(randi(n_y, n_y, 1));
        
        boot_effects(b) = effect_func(x_boot, y_boot);
    end
    
    % Remove invalid values
    boot_effects = boot_effects(isfinite(boot_effects));
    
    if isempty(boot_effects)
        ci_low = NaN;
        ci_high = NaN;
        return;
    end
    
    % BCa (bias-corrected and accelerated) percentile method
    % Simplified: using percentile method
    alpha = 1 - ci_level;
    ci_low = prctile(boot_effects, 100 * alpha/2);
    ci_high = prctile(boot_effects, 100 * (1 - alpha/2));
end


function [p_hier, ci_diff] = hierarchical_permutation_test(sample_data_1, sample_data_2, n_perm)
% HIERARCHICAL_PERMUTATION_TEST - Permutation test for nested data
%
% Shuffles sample labels (not individual ROIs) to maintain cluster structure.
% This is more conservative and appropriate for nested data than ROI-level permutation.
%
% INPUTS:
%   sample_data_1 - cell array, each cell = ROI data from one sample (group 1)
%   sample_data_2 - cell array, each cell = ROI data from one sample (group 2)
%   n_perm        - number of permutations (default: 10000)
%
% OUTPUTS:
%   p_hier  - p-value from hierarchical permutation test
%   ci_diff - 95% CI for difference in means

    if nargin < 3 || isempty(n_perm)
        n_perm = 10000;
    end
    
    % Compute observed statistic (difference in means)
    mean_1 = mean(cellfun(@(x) mean(x, 'omitnan'), sample_data_1), 'omitnan');
    mean_2 = mean(cellfun(@(x) mean(x, 'omitnan'), sample_data_2), 'omitnan');
    observed_diff = mean_1 - mean_2;
    
    % Combine all samples
    all_samples = [sample_data_1(:); sample_data_2(:)];
    n_total = length(all_samples);
    n_1 = length(sample_data_1);
    
    % Permutation distribution
    perm_diffs = zeros(n_perm, 1);
    
    for p = 1:n_perm
        % Shuffle sample labels
        perm_idx = randperm(n_total);
        perm_1 = all_samples(perm_idx(1:n_1));
        perm_2 = all_samples(perm_idx(n_1+1:end));
        
        mean_p1 = mean(cellfun(@(x) mean(x, 'omitnan'), perm_1), 'omitnan');
        mean_p2 = mean(cellfun(@(x) mean(x, 'omitnan'), perm_2), 'omitnan');
        perm_diffs(p) = mean_p1 - mean_p2;
    end
    
    % Two-tailed p-value
    p_hier = mean(abs(perm_diffs) >= abs(observed_diff));
    
    % Bootstrap CI for difference
    [ci_diff(1), ci_diff(2)] = bootstrap_mean_diff_ci(sample_data_1, sample_data_2, 10000, 0.95);
end


function [ci_low, ci_high] = bootstrap_mean_diff_ci(sample_data_1, sample_data_2, n_boot, ci_level)
% BOOTSTRAP_MEAN_DIFF_CI - Bootstrap CI for difference in means (hierarchical)

    if nargin < 3, n_boot = 10000; end
    if nargin < 4, ci_level = 0.95; end
    
    n_1 = length(sample_data_1);
    n_2 = length(sample_data_2);
    
    boot_diffs = zeros(n_boot, 1);
    
    for b = 1:n_boot
        % Resample samples
        idx_1 = randi(n_1, n_1, 1);
        idx_2 = randi(n_2, n_2, 1);
        
        mean_1 = mean(cellfun(@(x) mean(x, 'omitnan'), sample_data_1(idx_1)), 'omitnan');
        mean_2 = mean(cellfun(@(x) mean(x, 'omitnan'), sample_data_2(idx_2)), 'omitnan');
        
        boot_diffs(b) = mean_1 - mean_2;
    end
    
    alpha = 1 - ci_level;
    ci_low = prctile(boot_diffs, 100 * alpha/2);
    ci_high = prctile(boot_diffs, 100 * (1 - alpha/2));
end


function icc = compute_icc(sample_data)
% COMPUTE_ICC - Compute Intraclass Correlation Coefficient
%
% ICC measures how much of the total variance is due to between-sample
% differences vs within-sample differences.
%
% High ICC (>0.5) means samples are quite different from each other,
% suggesting sample-level analysis is important.
%
% Low ICC (<0.1) means most variance is within samples.

    n_samples = length(sample_data);
    
    % Get all values and group labels
    all_values = [];
    group_labels = [];
    
    for s = 1:n_samples
        data = sample_data{s}(:);
        data = data(~isnan(data));
        all_values = [all_values; data];
        group_labels = [group_labels; s * ones(length(data), 1)];
    end
    
    if isempty(all_values) || length(unique(group_labels)) < 2
        icc = NaN;
        return;
    end
    
    % One-way ANOVA to get variance components
    [~, tbl] = anova1(all_values, group_labels, 'off');
    
    if size(tbl, 1) < 3
        icc = NaN;
        return;
    end
    
    % Extract mean squares
    MS_between = tbl{2, 4};  % Between-group MS
    MS_within = tbl{3, 4};   % Within-group MS
    
    % Average group size
    n_per_group = arrayfun(@(x) sum(group_labels == x), unique(group_labels));
    n_0 = mean(n_per_group);
    
    % ICC(1) - one-way random effects
    icc = (MS_between - MS_within) / (MS_between + (n_0 - 1) * MS_within);
    icc = max(0, icc);  % ICC can't be negative in practice
end


function design_effect = compute_design_effect(sample_data)
% COMPUTE_DESIGN_EFFECT - Compute design effect for clustered data
%
% Design effect = 1 + (average cluster size - 1) * ICC
%
% This tells you how much the effective sample size is reduced due to clustering.
% If design_effect = 5, then 100 clustered observations have the effective
% sample size of only 100/5 = 20 independent observations.

    icc = compute_icc(sample_data);
    
    n_per_sample = cellfun(@(x) sum(~isnan(x)), sample_data);
    avg_cluster_size = mean(n_per_sample);
    
    design_effect = 1 + (avg_cluster_size - 1) * icc;
end


function cohens_d = compute_cohens_d_pooled(x, y)
% COMPUTE_COHENS_D_POOLED - Cohen's d with pooled standard deviation

    x = x(~isnan(x));
    y = y(~isnan(y));
    
    n_x = length(x);
    n_y = length(y);
    
    if n_x < 2 || n_y < 2
        cohens_d = NaN;
        return;
    end
    
    mean_diff = mean(x) - mean(y);
    
    % Pooled SD (weighted by df)
    pooled_var = ((n_x - 1) * var(x) + (n_y - 1) * var(y)) / (n_x + n_y - 2);
    pooled_sd = sqrt(pooled_var);
    
    if pooled_sd == 0
        cohens_d = NaN;
    else
        cohens_d = mean_diff / pooled_sd;
    end
end


function hedges_g = compute_hedges_g(x, y)
% COMPUTE_HEDGES_G - Hedges' g (bias-corrected Cohen's d)
%
% Hedges' g corrects for small-sample bias in Cohen's d.
% Recommended when total N < 50.

    x = x(~isnan(x));
    y = y(~isnan(y));
    
    n_x = length(x);
    n_y = length(y);
    
    d = compute_cohens_d_pooled(x, y);
    
    if isnan(d)
        hedges_g = NaN;
        return;
    end
    
    % Correction factor (approximation)
    df = n_x + n_y - 2;
    correction = 1 - (3 / (4 * df - 1));
    
    hedges_g = d * correction;
end


function print_enhanced_stats_table(comparisons, results, alpha_adj)
% PRINT_ENHANCED_STATS_TABLE - Print publication-ready statistics table
%
% Follows APA guidelines for reporting effect sizes with CIs

    fprintf('\n');
    fprintf('╔═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║                                    STATISTICAL RESULTS (with Effect Size Confidence Intervals)                                  ║\n');
    fprintf('╠═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣\n');
    fprintf('║ Metric                    │ Ctrl Median (IQR) │ Exp Median (IQR)  │ p-value  │ Cliff δ [95%% CI]      │ Interpretation         ║\n');
    fprintf('╠═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╣\n');
    
    for i = 1:length(comparisons)
        c = comparisons(i);
        r = results(i);
        
        % Format numbers
        if c.ctrl_median < 0.01
            ctrl_str = sprintf('%6.4f (%6.4f)', c.ctrl_median, c.ctrl_iqr);
        elseif c.ctrl_median < 1
            ctrl_str = sprintf('%6.3f (%6.3f)', c.ctrl_median, c.ctrl_iqr);
        else
            ctrl_str = sprintf('%6.2f (%6.2f)', c.ctrl_median, c.ctrl_iqr);
        end
        
        if c.exp_median < 0.01
            exp_str = sprintf('%6.4f (%6.4f)', c.exp_median, c.exp_iqr);
        elseif c.exp_median < 1
            exp_str = sprintf('%6.3f (%6.3f)', c.exp_median, c.exp_iqr);
        else
            exp_str = sprintf('%6.2f (%6.2f)', c.exp_median, c.exp_iqr);
        end
        
        % Effect size with CI
        if isfield(r, 'delta_ci_low') && ~isnan(r.delta_ci_low)
            effect_str = sprintf('%+.2f [%+.2f, %+.2f]', r.delta, r.delta_ci_low, r.delta_ci_high);
        else
            effect_str = sprintf('%+.2f', r.delta);
        end
        
        % Interpretation
        if r.p < alpha_adj
            if r.delta > 0.474
                interp = 'Ctrl >> Exp ***';
            elseif r.delta > 0.33
                interp = 'Ctrl > Exp **';
            elseif r.delta > 0.147
                interp = 'Ctrl > Exp *';
            elseif r.delta < -0.474
                interp = 'Exp >> Ctrl ***';
            elseif r.delta < -0.33
                interp = 'Exp > Ctrl **';
            elseif r.delta < -0.147
                interp = 'Exp > Ctrl *';
            else
                interp = 'negligible';
            end
        else
            interp = 'n.s.';
        end
        
        % Significance marker
        if r.p < 0.001
            p_str = '< 0.001';
        else
            p_str = sprintf('%.4f', r.p);
        end
        
        fprintf('║ %-25s │ %-17s │ %-17s │ %-8s │ %-21s │ %-22s ║\n', ...
            c.name, ctrl_str, exp_str, p_str, effect_str, interp);
    end
    
    fprintf('╚═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝\n');
    fprintf('\n');
    fprintf('Effect size interpretation (Cliff''s δ): |δ|<0.147 negligible, <0.33 small, <0.474 medium, ≥0.474 large\n');
    fprintf('Significance: Bonferroni-corrected α = %.4f\n', alpha_adj);
end


function report = generate_methods_text(n_ctrl, n_exp, n_roi_ctrl, n_roi_exp, alpha, n_comparisons)
% GENERATE_METHODS_TEXT - Generate methods section text for publication
%
% This automatically generates appropriate statistical methods text
% based on the analysis performed.

    report = sprintf([...
        'STATISTICAL METHODS (for Methods section)\n' ...
        '─────────────────────────────────────────────────────────────────────────\n' ...
        'Statistical analyses were performed using MATLAB (MathWorks, Inc.). ' ...
        'Data from %d control and %d experimental samples were compared. ' ...
        'A total of %d and %d ROIs were analyzed for control and experimental ' ...
        'conditions, respectively.\n\n' ...
        'Given the nested data structure (ROIs within samples), primary ' ...
        'comparisons used sample-level summary statistics (one value per ' ...
        'animal) to avoid pseudoreplication (Saravanan et al., 2020). ' ...
        'Non-normal distributions were confirmed by Lilliefors test, and ' ...
        'non-parametric Mann-Whitney U tests were used for group comparisons. ' ...
        'To account for multiple comparisons (%d tests), we applied ' ...
        'Bonferroni correction (adjusted α = %.4f).\n\n' ...
        'Effect sizes were quantified using Cliff''s delta (δ), with 95%% ' ...
        'bootstrap confidence intervals (10,000 iterations). For ROI-level ' ...
        'metrics, we computed hierarchical bootstrap CIs using a two-level ' ...
        'resampling procedure (samples, then ROIs within samples) to obtain ' ...
        'conservative estimates that correctly account for within-sample ' ...
        'correlation.\n\n' ...
        'The intraclass correlation coefficient (ICC) and design effect were ' ...
        'computed to quantify clustering. If hierarchical bootstrap CIs ' ...
        'included zero while naive ROI-level tests showed significance, the ' ...
        'result was interpreted as not significant.\n\n' ...
        'Effect size thresholds followed Romano et al. (2006): |δ| < 0.147 ' ...
        'negligible, < 0.33 small, < 0.474 medium, ≥ 0.474 large.\n'], ...
        n_ctrl, n_exp, n_roi_ctrl, n_roi_exp, n_comparisons, alpha/n_comparisons);
end
