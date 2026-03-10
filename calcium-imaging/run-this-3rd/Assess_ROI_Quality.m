function Assess_ROI_Quality(varargin)
%% ASSESS_ROI_QUALITY v2.0 - Post-Pipeline Signal Quality Assessment
% =========================================================================
% Run this AFTER the main pipeline to assess signal quality of all samples
%
% USAGE:
%   Assess_ROI_Quality()           % Uses defaults + folder picker
%   Assess_ROI_Quality(cfg)        % Uses custom config struct
%
%   cfg = assess_quality_config();
%   cfg.conditions(1).name = 'Control';
%   cfg.conditions(1).roi_dir = '/path/to/roi_analysis';
%   cfg.conditions(1).color = [0.2 0.4 0.8];
%   cfg.output_dir = '/path/to/output';
%   Assess_ROI_Quality(cfg);
%
% v2.0 UPDATE (February 2026):
%   - Now reads quality_metrics from v4.5 pipeline output
%   - Uses pairwise SYNCHRONY (correct metric) instead of ROI-mean correlation
%   - Assesses dFF_raw (pre-regression) for contamination severity
%   - Recalibrated thresholds for synchrony metric
%   - Falls back to computing metrics for v4.4 data
%   - Removed editorial conclusions from auto-generated README
%
% PURPOSE:
%   The main pipeline (COMPLETE_PIPELINE_RUN_ME) handles ROI detection.
%   This script evaluates the QUALITY of detected signals to:
%   1. Flag samples with severe artifacts (extreme dF/F values)
%   2. Flag samples with high neuropil contamination (pairwise synchrony)
%   3. Create eligibility lists for downstream analysis
%   4. Generate publication-ready summary comparing conditions
%
% QUALITY CATEGORIES (based on PAIRWISE SYNCHRONY):
%   ARTIFACT       - Extreme dF/F values (>5 or <-2), sample unusable
%   SEVERE_NEUROPIL- Very high synchrony (>0.5), sample unusable
%   HIGH_NEUROPIL  - High synchrony (>0.3), needs global regression
%   MILD_NEUROPIL  - Moderate synchrony (0.15-0.3), acceptable
%   OK             - Clean signal (<0.15 pairwise synchrony)
%
% NOTE: Pairwise synchrony values are typically LOWER than ROI-mean correlation.
%   Synchrony 0.15 ~ ROI-mean corr 0.4
%   Synchrony 0.30 ~ ROI-mean corr 0.6
%   Synchrony 0.50 ~ ROI-mean corr 0.7+
%
% OUTPUTS (saved to output directory):
%   - sample_quality_assessment.csv    : All samples with detailed metrics
%   - eligible_samples_control.txt     : List of usable Control samples
%   - eligible_samples_experimental.txt: List of usable Experimental samples
%   - quality_comparison_summary.csv   : Summary statistics
%   - quality_diagnostic_figures/      : Per-sample diagnostic plots
%   - README_quality_assessment.txt    : Documentation
%
% CONFIGURATION:
%   All parameters are set via the config struct. See
%   assess_quality_config() for the full list and documentation.
%
% See also: assess_quality_config, COMPLETE_PIPELINE_RUN_ME_v4_5
% =========================================================================

%% ========================================================================
%  CONFIGURATION
% =========================================================================

if nargin >= 1 && isstruct(varargin{1})
    cfg = merge_with_defaults(varargin{1});
elseif nargin >= 1
    error(['Unexpected argument (type: %s). Pass a config struct from\n' ...
           'assess_quality_config(), or call with no arguments for folder picker.'], ...
           class(varargin{1}));
else
    cfg = assess_quality_config();
end

% If no conditions defined, prompt for folders
if isempty(cfg.conditions) || ~isfield(cfg.conditions, 'name') || isempty(cfg.conditions(1).name)
    cfg.conditions = prompt_for_conditions();
end

% If no output dir, prompt
if isempty(cfg.output_dir)
    if usejava('desktop')
        cfg.output_dir = uigetdir('', 'Select output directory for quality assessment');
        if cfg.output_dir == 0
            error('No output directory selected. Exiting.');
        end
    else
        error(['No output_dir specified. Set cfg.output_dir before calling.\n' ...
               'See also: assess_quality_config']);
    end
end

CONDITIONS = cfg.conditions;
OUTPUT_DIR = cfg.output_dir;
THRESHOLDS = cfg.thresholds;
OPTIONS = cfg.options;

%% ========================================================================
%  INITIALIZATION
% =========================================================================

fprintf('\n');
fprintf('==================================================================\n');
fprintf('  POST-PIPELINE ROI QUALITY ASSESSMENT v2.0\n');
fprintf('  Using Pairwise Synchrony (Pipeline v4.5 Compatible)\n');
fprintf('==================================================================\n');
fprintf('\n');
fprintf('Date: %s\n', datetime('now'));
fprintf('\n');

fprintf('=== QUALITY THRESHOLDS (Pairwise Synchrony) ===\n');
fprintf('  Artifact: dF/F > %.1f or < %.1f\n', THRESHOLDS.max_dff, THRESHOLDS.min_dff);
fprintf('  Severe neuropil: synchrony > %.2f\n', THRESHOLDS.severe_neuropil);
fprintf('  High neuropil: synchrony > %.2f\n', THRESHOLDS.high_neuropil);
fprintf('  Mild neuropil: synchrony > %.2f\n', THRESHOLDS.mild_neuropil);
fprintf('  Minimum ROIs: %d\n', THRESHOLDS.min_rois);
fprintf('\n');

% Create output directory
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('Created output directory: %s\n', OUTPUT_DIR);
end

if OPTIONS.generate_diagnostic_figures
    fig_dir = fullfile(OUTPUT_DIR, 'quality_diagnostic_figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
end

%% ========================================================================
%  ASSESS ALL SAMPLES
% =========================================================================

all_samples = table();

for c = 1:length(CONDITIONS)
    cond = CONDITIONS(c);
    
    fprintf('\n');
    fprintf('==================================================================\n');
    fprintf('  CONDITION: %s\n', cond.name);
    fprintf('==================================================================\n');
    
    if ~exist(cond.roi_dir, 'dir')
        fprintf('  [SKIP] Directory not found: %s\n', cond.roi_dir);
        continue;
    end
    
    mat_files = dir(fullfile(cond.roi_dir, '*_rois.mat'));
    
    if isempty(mat_files)
        fprintf('  [SKIP] No ROI files found\n');
        continue;
    end
    
    fprintf('  Found %d ROI files\n\n', length(mat_files));
    
    n_diagnostic_figs = 0;
    
    for f = 1:length(mat_files)
        mat_path = fullfile(mat_files(f).folder, mat_files(f).name);
        sample_name = strrep(mat_files(f).name, '_rois.mat', '');
        
        fprintf('  --- %s ---\n', sample_name);
        
        % Load data
        try
            loaded = load(mat_path);
            if isfield(loaded, 'roi_data')
                roi_data = loaded.roi_data;
            else
                fnames = fieldnames(loaded);
                roi_data = loaded.(fnames{1});
            end
        catch ME
            fprintf('    [ERROR] Failed to load: %s\n', ME.message);
            continue;
        end
        
        % Initialize row
        row = table();
        row.sample = {sample_name};
        row.condition = {cond.name};
        
        % Get ROI count
        if isfield(roi_data, 'roi_centers')
            n_rois = size(roi_data.roi_centers, 1);
        else
            n_rois = 0;
        end
        row.n_rois = n_rois;
        
        % Handle no ROIs
        if n_rois == 0
            row.max_dff = NaN;
            row.min_dff = NaN;
            row.synchrony = NaN;
            row.synchrony_after = NaN;
            row.neuropil_severity = {'NONE'};
            row.regression_applied = false;
            row.mean_snr = NaN;
            row.mean_event_rate = NaN;
            row.status = {'NO_ROIS'};
            row.eligible = false;
            row.reason = {'No ROIs detected'};
            row.needs_regression = false;
            row.pipeline_version = {'unknown'};
            
            all_samples = [all_samples; row];
            fprintf('    Status: NO_ROIS\n');
            continue;
        end
        
        % =================================================================
        % CHECK FOR v4.5 QUALITY METRICS (preferred)
        % =================================================================
        
        artifact_flag = false;
        
        if isfield(roi_data, 'quality_metrics')
            % v4.5 data - use stored metrics
            qm = roi_data.quality_metrics;
            row.pipeline_version = {'v4.5'};
            
            % Get synchrony (the correct metric)
            if isfield(qm, 'synchrony_before')
                row.synchrony = qm.synchrony_before;
            else
                row.synchrony = NaN;
            end
            
            if isfield(qm, 'synchrony_after')
                row.synchrony_after = qm.synchrony_after;
            else
                row.synchrony_after = NaN;
            end
            
            % Get artifact info
            if isfield(qm, 'max_dff')
                row.max_dff = qm.max_dff;
            else
                row.max_dff = NaN;
            end
            
            if isfield(qm, 'min_dff')
                row.min_dff = qm.min_dff;
            else
                row.min_dff = NaN;
            end
            
            % Get neuropil severity from pipeline
            if isfield(qm, 'neuropil_severity')
                row.neuropil_severity = {qm.neuropil_severity};
            else
                row.neuropil_severity = {'UNKNOWN'};
            end
            
            % Check if regression was applied
            if isfield(qm, 'global_regression_applied')
                row.regression_applied = qm.global_regression_applied;
            else
                row.regression_applied = false;
            end
            
            % Check artifact flag
            if isfield(qm, 'artifact_flag')
                artifact_flag = qm.artifact_flag;
            end
            
            fprintf('    [v4.5] synchrony=%.3f, regression=%s\n', ...
                row.synchrony, mat2str(row.regression_applied));
            
        else
            % =============================================================
            % FALLBACK FOR v4.4 DATA: Compute metrics manually
            % =============================================================
            
            row.pipeline_version = {'v4.4'};
            row.regression_applied = false;
            row.synchrony_after = NaN;
            
            fprintf('    [v4.4 fallback] Computing metrics...\n');
            
            % Get traces - use dFF_raw if available (pre-regression)
            if isfield(roi_data, 'dFF_raw')
                dFF = double(roi_data.dFF_raw);
            else
                dFF = double(roi_data.dFF);
            end
            
            % Pipeline ALWAYS stores dFF as [nROIs × nFrames]
            % Transpose to [nFrames × nROIs] for corrcoef (each column = 1 ROI)
            dFF = dFF';
            
            % Calculate dF/F extremes
            row.max_dff = max(dFF(:));
            row.min_dff = min(dFF(:));
            
            % Check for artifacts
            if row.max_dff > THRESHOLDS.max_dff || row.min_dff < THRESHOLDS.min_dff
                artifact_flag = true;
            end
            
            % Calculate PAIRWISE SYNCHRONY (correct metric)
            if size(dFF, 2) >= 2
                corr_matrix = corrcoef(dFF);
                upper_tri = triu(true(size(corr_matrix)), 1);
                pairwise_corrs = corr_matrix(upper_tri);
                row.synchrony = mean(abs(pairwise_corrs), 'omitnan');
            else
                row.synchrony = NaN;
            end
            
            % Classify neuropil severity
            if row.synchrony > THRESHOLDS.severe_neuropil
                row.neuropil_severity = {'SEVERE'};
            elseif row.synchrony > THRESHOLDS.high_neuropil
                row.neuropil_severity = {'HIGH'};
            elseif row.synchrony > THRESHOLDS.mild_neuropil
                row.neuropil_severity = {'MILD'};
            else
                row.neuropil_severity = {'LOW'};
            end
        end
        
        % Get other metrics
        if isfield(roi_data, 'snr')
            row.mean_snr = mean(roi_data.snr);
        else
            row.mean_snr = NaN;
        end
        
        if isfield(roi_data, 'event_rate')
            row.mean_event_rate = mean(roi_data.event_rate);
        else
            row.mean_event_rate = NaN;
        end
        
        % -----------------------------------------------------------------
        % DETERMINE STATUS
        % -----------------------------------------------------------------
        
        if artifact_flag || row.max_dff > THRESHOLDS.max_dff || row.min_dff < THRESHOLDS.min_dff
            row.status = {'ARTIFACT'};
            row.eligible = false;
            row.reason = {sprintf('Extreme dF/F (max=%.1f, min=%.1f)', row.max_dff, row.min_dff)};
            row.needs_regression = false;
            
        elseif n_rois < THRESHOLDS.min_rois
            row.status = {'TOO_FEW_ROIS'};
            row.eligible = false;
            row.reason = {sprintf('Only %d ROIs (need >= %d)', n_rois, THRESHOLDS.min_rois)};
            row.needs_regression = false;
            
        elseif row.synchrony > THRESHOLDS.severe_neuropil
            row.status = {'SEVERE_NEUROPIL'};
            row.eligible = false;
            row.reason = {sprintf('Severe neuropil (synchrony=%.3f)', row.synchrony)};
            row.needs_regression = true;
            
        elseif row.synchrony > THRESHOLDS.high_neuropil
            row.status = {'HIGH_NEUROPIL'};
            row.eligible = true;
            row.reason = {sprintf('High neuropil (synchrony=%.3f)', row.synchrony)};
            row.needs_regression = true;
            
        elseif row.synchrony > THRESHOLDS.mild_neuropil
            row.status = {'MILD_NEUROPIL'};
            row.eligible = true;
            row.reason = {sprintf('Mild neuropil (synchrony=%.3f)', row.synchrony)};
            row.needs_regression = true;
            
        else
            row.status = {'OK'};
            row.eligible = true;
            row.reason = {'Clean signal'};
            row.needs_regression = false;
        end
        
        all_samples = [all_samples; row];
        
        % Print status
        if row.eligible
            status_sym = '+';
        else
            status_sym = 'X';
        end
        fprintf('    ROIs: %d | Synchrony: %.3f | Status: %s %s\n', ...
            n_rois, row.synchrony, status_sym, row.status{1});
        
        % Generate diagnostic figure
        if OPTIONS.generate_diagnostic_figures && n_diagnostic_figs < OPTIONS.max_diagnostic_figures
            if row.synchrony > THRESHOLDS.mild_neuropil || ~row.eligible
                if isfield(roi_data, 'dFF_raw')
                    dFF_viz = double(roi_data.dFF_raw);
                else
                    dFF_viz = double(roi_data.dFF);
                end
                
                % Pipeline stores [nROIs × nFrames], transpose for visualization
                dFF_viz = dFF_viz';
                
                generate_diagnostic_figure_v2(dFF_viz, fig_dir, sample_name, ...
                    cond.name, row, THRESHOLDS);
                n_diagnostic_figs = n_diagnostic_figs + 1;
            end
        end
    end
end

%% ========================================================================
%  GENERATE OUTPUTS
% =========================================================================

if isempty(all_samples)
    fprintf('\n[ERROR] No samples found to assess!\n');
    return;
end

fprintf('\n');
fprintf('==================================================================\n');
fprintf('  GENERATING OUTPUTS\n');
fprintf('==================================================================\n\n');

% 1. Full assessment CSV
assessment_path = fullfile(OUTPUT_DIR, 'sample_quality_assessment.csv');
writetable(all_samples, assessment_path);
fprintf('+ Saved: %s\n', assessment_path);

% 2. Eligible sample lists
for c = 1:length(CONDITIONS)
    cond_name = CONDITIONS(c).name;
    cond_rows = strcmp(all_samples.condition, cond_name);
    eligible_rows = cond_rows & all_samples.eligible;
    
    eligible_samples = all_samples.sample(eligible_rows);
    
    list_path = fullfile(OUTPUT_DIR, sprintf('eligible_samples_%s.txt', lower(strrep(cond_name, ' ', '_'))));
    fid = fopen(list_path, 'w');
    for i = 1:length(eligible_samples)
        fprintf(fid, '%s\n', eligible_samples{i});
    end
    fclose(fid);
    fprintf('+ Saved: %s (%d samples)\n', list_path, length(eligible_samples));
end

% 3. Quality comparison summary
summary = table();

for c = 1:length(CONDITIONS)
    cond_name = CONDITIONS(c).name;
    cond_rows = strcmp(all_samples.condition, cond_name);
    
    row = table();
    row.condition = {cond_name};
    row.total_samples = sum(cond_rows);
    row.eligible_samples = sum(cond_rows & all_samples.eligible);
    row.pct_eligible = 100 * row.eligible_samples / row.total_samples;
    
    row.n_ok = sum(cond_rows & strcmp(all_samples.status, 'OK'));
    row.n_mild_neuropil = sum(cond_rows & strcmp(all_samples.status, 'MILD_NEUROPIL'));
    row.n_high_neuropil = sum(cond_rows & strcmp(all_samples.status, 'HIGH_NEUROPIL'));
    row.n_severe_neuropil = sum(cond_rows & strcmp(all_samples.status, 'SEVERE_NEUROPIL'));
    row.n_artifact = sum(cond_rows & strcmp(all_samples.status, 'ARTIFACT'));
    row.n_no_rois = sum(cond_rows & strcmp(all_samples.status, 'NO_ROIS'));
    row.n_too_few_rois = sum(cond_rows & strcmp(all_samples.status, 'TOO_FEW_ROIS'));
    
    row.pct_artifact = 100 * row.n_artifact / row.total_samples;
    row.pct_clean = 100 * row.n_ok / row.total_samples;
    
    eligible_idx = cond_rows & all_samples.eligible;
    if any(eligible_idx)
        row.mean_synchrony = mean(all_samples.synchrony(eligible_idx), 'omitnan');
        row.mean_snr = mean(all_samples.mean_snr(eligible_idx), 'omitnan');
        row.mean_rois = mean(all_samples.n_rois(eligible_idx));
    else
        row.mean_synchrony = NaN;
        row.mean_snr = NaN;
        row.mean_rois = NaN;
    end
    
    summary = [summary; row];
end

summary_path = fullfile(OUTPUT_DIR, 'quality_comparison_summary.csv');
writetable(summary, summary_path);
fprintf('+ Saved: %s\n', summary_path);

% 4. README (objective, no editorial conclusions)
readme_path = fullfile(OUTPUT_DIR, 'README_quality_assessment.txt');
fid = fopen(readme_path, 'w');

fprintf(fid, '===============================================================================\n');
fprintf(fid, 'ROI QUALITY ASSESSMENT v2.0 - README\n');
fprintf(fid, '===============================================================================\n\n');

fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));

fprintf(fid, 'METRIC: PAIRWISE SYNCHRONY\n');
fprintf(fid, '-------------------------------------------------------------------------------\n');
fprintf(fid, 'This script uses PAIRWISE SYNCHRONY (mean absolute pairwise correlation)\n');
fprintf(fid, 'rather than ROI-mean correlation. This avoids circularity.\n\n');
fprintf(fid, 'Typical values:\n');
fprintf(fid, '  - Independent neurons: synchrony < 0.15\n');
fprintf(fid, '  - Mild shared signal: synchrony 0.15 - 0.30\n');
fprintf(fid, '  - Strong contamination: synchrony > 0.30\n');
fprintf(fid, '  - Severe (likely artifact): synchrony > 0.50\n\n');

fprintf(fid, 'THRESHOLDS\n');
fprintf(fid, '-------------------------------------------------------------------------------\n');
fprintf(fid, '  max_dff = %.1f\n', THRESHOLDS.max_dff);
fprintf(fid, '  min_dff = %.1f\n', THRESHOLDS.min_dff);
fprintf(fid, '  severe_neuropil = %.2f\n', THRESHOLDS.severe_neuropil);
fprintf(fid, '  high_neuropil = %.2f\n', THRESHOLDS.high_neuropil);
fprintf(fid, '  mild_neuropil = %.2f\n', THRESHOLDS.mild_neuropil);
fprintf(fid, '  min_rois = %d\n\n', THRESHOLDS.min_rois);

fprintf(fid, 'RESULTS\n');
fprintf(fid, '-------------------------------------------------------------------------------\n');

for c = 1:height(summary)
    fprintf(fid, '\n%s:\n', summary.condition{c});
    fprintf(fid, '  Total: %d | Eligible: %d (%.1f%%)\n', ...
        summary.total_samples(c), summary.eligible_samples(c), summary.pct_eligible(c));
    fprintf(fid, '  OK: %d | Artifact: %d | Severe neuropil: %d\n', ...
        summary.n_ok(c), summary.n_artifact(c), summary.n_severe_neuropil(c));
    fprintf(fid, '  Mean synchrony (eligible): %.3f\n', summary.mean_synchrony(c));
end

fprintf(fid, '\n');
fprintf(fid, 'COMPARISON (numbers only - interpret in your paper)\n');
fprintf(fid, '-------------------------------------------------------------------------------\n');

if height(summary) >= 2
    ctrl_idx = strcmp(summary.condition, 'Control');
    exp_idx = strcmp(summary.condition, 'Experimental');
    
    if any(ctrl_idx) && any(exp_idx)
        fprintf(fid, '                        Control     Experimental\n');
        fprintf(fid, '  Eligible (%%)          %5.1f       %5.1f\n', ...
            summary.pct_eligible(ctrl_idx), summary.pct_eligible(exp_idx));
        fprintf(fid, '  Artifact (%%)          %5.1f       %5.1f\n', ...
            summary.pct_artifact(ctrl_idx), summary.pct_artifact(exp_idx));
        fprintf(fid, '  Clean (%%)             %5.1f       %5.1f\n', ...
            summary.pct_clean(ctrl_idx), summary.pct_clean(exp_idx));
        fprintf(fid, '  Mean synchrony        %5.3f       %5.3f\n', ...
            summary.mean_synchrony(ctrl_idx), summary.mean_synchrony(exp_idx));
    end
end

fprintf(fid, '\n===============================================================================\n');
fclose(fid);
fprintf('+ Saved: %s\n', readme_path);

%% ========================================================================
%  PRINT SUMMARY
% =========================================================================

fprintf('\n');
fprintf('==================================================================\n');
fprintf('  QUALITY ASSESSMENT COMPLETE\n');
fprintf('==================================================================\n');
fprintf('\n');

fprintf('%-15s %8s %8s %8s %8s %10s\n', 'Condition', 'Total', 'Eligible', 'Artifact', 'Clean', 'Synchrony');
fprintf('%s\n', repmat('-', 1, 65));

for c = 1:height(summary)
    fprintf('%-15s %8d %8d %8d %8d %10.3f\n', ...
        summary.condition{c}, ...
        summary.total_samples(c), ...
        summary.eligible_samples(c), ...
        summary.n_artifact(c), ...
        summary.n_ok(c), ...
        summary.mean_synchrony(c));
end

fprintf('\nOutput: %s\n\n', OUTPUT_DIR);

end  % end of main function Assess_ROI_Quality

%% ========================================================================
%  CONFIG HELPERS
% =========================================================================

function cfg = merge_with_defaults(user_cfg)
cfg = assess_quality_config();
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
                end
            end
        else
            cfg.(user_fields{i}) = user_cfg.(user_fields{i});
        end
    end
end
end

function conditions = prompt_for_conditions()
if ~usejava('desktop')
    error(['No conditions defined and no GUI available.\n' ...
           'See also: assess_quality_config']);
end
default_colors = {[0.2 0.4 0.8], [0.9 0.5 0.1], [0.2 0.7 0.3], [0.6 0.2 0.7]};
conditions = struct('name', {}, 'roi_dir', {}, 'color', {});
condIdx = 0;
while true
    condIdx = condIdx + 1;
    answer = questdlg(sprintf('Add condition %d?', condIdx), ...
        'Add Condition', 'Yes', 'No (done)', 'Yes');
    if ~strcmp(answer, 'Yes'), break; end
    name = inputdlg(sprintf('Condition %d name:', condIdx), 'Condition Name', 1, {sprintf('Condition_%d', condIdx)});
    if isempty(name), break; end
    roiDir = uigetdir('', sprintf('Select ROI analysis folder for %s', name{1}));
    if roiDir == 0, break; end
    conditions(condIdx).name = name{1};
    conditions(condIdx).roi_dir = roiDir;
    conditions(condIdx).color = default_colors{min(condIdx, length(default_colors))};
end
if isempty(conditions), error('No conditions defined.'); end
end

%% ========================================================================
%  HELPER: DIAGNOSTIC FIGURE
% =========================================================================

function generate_diagnostic_figure_v2(dFF, fig_dir, sample_name, cond_name, metrics, thresholds)

n_frames = size(dFF, 1);
n_rois = size(dFF, 2);

fig = figure('Position', [50 50 1400 600], 'Visible', 'off');

% Panel 1: Mean trace
subplot(2, 2, 1);
mean_trace = mean(dFF, 2);
plot(mean_trace, 'k', 'LineWidth', 1.5);
xlabel('Frame'); ylabel('dF/F');
title(sprintf('Mean trace (%d ROIs)', n_rois));
if metrics.synchrony > thresholds.high_neuropil
    set(gca, 'Color', [1 0.95 0.95]);
end

% Panel 2: Sample traces
subplot(2, 2, 2);
n_show = min(5, n_rois);
idx = randperm(n_rois, n_show);
hold on;
colors = lines(n_show);
for i = 1:n_show
    plot(dFF(:, idx(i)) + (i-1)*0.3, 'Color', colors(i,:), 'LineWidth', 0.5);
end
xlabel('Frame'); ylabel('dF/F (offset)');
title('Sample traces');

% Panel 3: Pairwise correlation histogram
subplot(2, 2, 3);
if n_rois >= 2
    corr_matrix = corrcoef(dFF);
    upper_tri = triu(true(n_rois), 1);
    pairwise_corrs = corr_matrix(upper_tri);
    
    histogram(pairwise_corrs, 30, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    hold on;
    xline(thresholds.mild_neuropil, 'g--', 'LineWidth', 2);
    xline(thresholds.high_neuropil, 'y--', 'LineWidth', 2);
    xline(thresholds.severe_neuropil, 'r--', 'LineWidth', 2);
    xline(metrics.synchrony, 'k-', 'LineWidth', 2);
    xlabel('Pairwise correlation'); ylabel('Count');
    title(sprintf('Pairwise synchrony: %.3f', metrics.synchrony));
    xlim([-0.5, 1]);
end

% Panel 4: Summary
subplot(2, 2, 4);
axis off;

summary_text = {
    sprintf('%s - %s', cond_name, sample_name), ...
    sprintf('Pipeline: %s', metrics.pipeline_version{1}), ...
    '', ...
    sprintf('Status: %s', metrics.status{1}), ...
    sprintf('Eligible: %s', mat2str(metrics.eligible)), ...
    '', ...
    sprintf('ROIs: %d', metrics.n_rois), ...
    sprintf('Synchrony: %.3f', metrics.synchrony), ...
    sprintf('Max dF/F: %.2f', metrics.max_dff), ...
    sprintf('Mean SNR: %.2f', metrics.mean_snr), ...
    '', ...
    sprintf('Reason: %s', metrics.reason{1})
};

text(0.1, 0.9, summary_text, 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'FontName', 'FixedWidth');

sgtitle(sprintf('%s - %s: Quality Diagnostic', cond_name, sample_name), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

saveas(fig, fullfile(fig_dir, sprintf('%s_%s_diagnostic.png', cond_name, sample_name)));
close(fig);

end