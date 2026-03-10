function ROI_REFINEMENT_v3_2(varargin)
%% ROI_REFINEMENT v3.2 - Strict Quality Control for Calcium Imaging
% =========================================================================
% COMPLETE VERSION - All visualization functions included
% UPDATED FOR PIPELINE v4.5
%
% USAGE:
%   ROI_REFINEMENT_v3_2()           % Uses defaults + folder picker
%   ROI_REFINEMENT_v3_2(cfg)        % Uses custom config struct
%
%   cfg = roi_refinement_config();
%   cfg.conditions(1).name = 'Control';
%   cfg.conditions(1).roi_analysis_dir = '/path/to/roi_analysis';
%   cfg.conditions(1).output_dir = '/path/to/refined';
%   ROI_REFINEMENT_v3_2(cfg);
%
% v3.2 UPDATE (February 2026):
%   - RECALIBRATED for wide-field imaging (was too aggressive)
%   - max_snr raised 50 → 150 (wide-field has higher baselines)
%   - Sustained high signal check relaxed: 5x→8x noise, 15%→30% threshold
%   - Kinetics thresholds widened for broad wide-field transients
%   - Shape criteria relaxed slightly (compactness 0.25→0.20, extent 0.40→0.35)
%
% v3.1 UPDATE (February 2026):
%   - Added dFF_raw passthrough for v4.5 compatibility
%   - Added quality_metrics passthrough (global_corr, neuropil_severity, etc.)
%
% v3.0 NEW CRITERIA (February 2026):
%
% SHAPE CRITERIA:
%   - Edge contact: Remove if >20% of ROI perimeter touches brain mask edge
%   - Compactness: Remove if shape is too irregular (tails, branches)
%   - Bounding box: Remove if ROI fills <40% or >95% of its bounding box
%   - Aspect ratio: Remove if length/width > 3 (too elongated)
%
% TRACE CRITERIA (PHYSIOLOGICAL VALIDATION):
%   - Negative transients: Remove if large dips below baseline
%   - Rise/decay kinetics: GCaMP8f should have fast rise, slow decay
%   - Step artifacts: Remove if instantaneous jumps (1-frame changes)
%   - Sustained plateaus: Remove if signal stays elevated without decay
%
% INHERITED FROM v2.3:
%   - Pairwise correlation > 0.7
%   - Spatial proximity < 5 px
%   - Solidity < 0.5
%   - Eccentricity > 0.95
%   - Event count = 0
%   - SNR > 50 (extreme artifacts)
%   - Baseline drift check
%
% OUTPUTS:
%   - Refined .mat files (same structure, fewer ROIs)
%   - ROI overlay image (all ROIs on brain with labels)
%   - Per-ROI gallery pages (30-second segmented traces)
%   - Refined analysis figures with adaptive heatmaps
%   - Rejection examples figure (shows WHY ROIs were removed)
%   - Refinement log (what was removed and why)
%   - Summary CSV ready for group comparison
%
% CONFIGURATION:
%   All parameters are set via the config struct. See
%   roi_refinement_config() for the full list and documentation.
%
% DEPENDENCIES:
%   - MATLAB R2020b or later
%   - Image Processing Toolbox (regionprops, bwperim)
%   - Signal Processing Toolbox (findpeaks)
%
% See also: roi_refinement_config, COMPLETE_PIPELINE_RUN_ME_v4_5
% =========================================================================

%% ════════════════════════════════════════════════════════════════════════
%  CONFIGURATION
% ════════════════════════════════════════════════════════════════════════

if nargin >= 1 && isstruct(varargin{1})
    cfg = merge_with_defaults(varargin{1});
elseif nargin >= 1
    error(['Unexpected argument (type: %s). Pass a config struct from\n' ...
           'roi_refinement_config(), or call with no arguments for folder picker.'], ...
           class(varargin{1}));
else
    cfg = roi_refinement_config();
end

% If no conditions defined, prompt for folders
if isempty(cfg.conditions) || ~isfield(cfg.conditions, 'name') || isempty(cfg.conditions(1).name)
    cfg.conditions = prompt_for_conditions();
end

CONDITIONS = cfg.conditions;

% Derive summary output dir if not set
if isempty(cfg.summary_output_dir)
    SUMMARY_OUTPUT_DIR = fullfile(fileparts(CONDITIONS(1).output_dir), 'ROI_Refinement_Summary');
else
    SUMMARY_OUTPUT_DIR = cfg.summary_output_dir;
end

% Imaging parameters
FRAME_RATE = cfg.frame_rate;
PIXEL_SIZE_UM = cfg.pixel_size_um;

%% ════════════════════════════════════════════════════════════════════════
%  REFINEMENT CRITERIA (from config)
% ════════════════════════════════════════════════════════════════════════

REFINE = struct();

% Shape criteria
REFINE.max_edge_contact_fraction = cfg.max_edge_contact_fraction;
REFINE.min_compactness = cfg.min_compactness;
REFINE.min_extent = cfg.min_extent;
REFINE.max_extent = cfg.max_extent;
REFINE.max_aspect_ratio = cfg.max_aspect_ratio;
REFINE.min_solidity = cfg.min_solidity;
REFINE.max_eccentricity = cfg.max_eccentricity;

% Trace / kinetics criteria
REFINE.max_negative_zscore = cfg.max_negative_zscore;
REFINE.max_rise_decay_ratio = cfg.max_rise_decay_ratio;
REFINE.max_single_frame_jump = cfg.max_single_frame_jump;
REFINE.max_plateau_duration_s = cfg.max_plateau_duration_s;
REFINE.plateau_threshold = cfg.plateau_threshold;
REFINE.min_decay_tau_s = cfg.min_decay_tau_s;
REFINE.max_decay_tau_s = cfg.max_decay_tau_s;

% Correlation / duplicate criteria
REFINE.max_pairwise_corr = cfg.max_pairwise_corr;
REFINE.corr_keep_higher_snr = cfg.corr_keep_higher_snr;
REFINE.min_center_distance_px = cfg.min_center_distance_px;

% Activity criteria
REFINE.min_event_count = cfg.min_event_count;
REFINE.max_snr = cfg.max_snr;

% Baseline criteria
REFINE.require_baseline_return = cfg.require_baseline_return;
REFINE.baseline_window_s = cfg.baseline_window_s;
REFINE.baseline_return_threshold = cfg.baseline_return_threshold;

% Output options
REFINE.generate_gallery = cfg.generate_gallery;
REFINE.generate_analysis_figures = cfg.generate_analysis_figures;
REFINE.generate_rejection_examples = cfg.generate_rejection_examples;
REFINE.gallery_segment_duration = cfg.gallery_segment_duration;
REFINE.skip_existing = cfg.skip_existing;
REFINE.organize_subfolders = cfg.organize_subfolders;

%% ════════════════════════════════════════════════════════════════════════
%  INITIALIZATION
% ════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║     ROI REFINEMENT v3.2 - RECALIBRATED FOR WIDE-FIELD                 ║\n');
fprintf('║     Shape + Kinetics Validation for Pipeline v4.5                ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Date: %s\n', datetime('now'));
fprintf('MATLAB: %s\n', version('-release'));
fprintf('Frame rate: %.2f Hz\n', FRAME_RATE);
fprintf('Pixel size: %.3f µm\n', PIXEL_SIZE_UM);
fprintf('\n');

fprintf('═══ v3.0 NEW SHAPE CRITERIA ═══\n');
fprintf('  Max edge contact: %.0f%% of perimeter\n', REFINE.max_edge_contact_fraction * 100);
fprintf('  Min compactness: %.2f (circle=1)\n', REFINE.min_compactness);
fprintf('  Extent range: %.2f - %.2f\n', REFINE.min_extent, REFINE.max_extent);
fprintf('  Max aspect ratio: %.1f\n', REFINE.max_aspect_ratio);
fprintf('  Min solidity: %.2f\n', REFINE.min_solidity);
fprintf('  Max eccentricity: %.2f\n', REFINE.max_eccentricity);
fprintf('\n');

fprintf('═══ v3.0 NEW KINETICS CRITERIA ═══\n');
fprintf('  Max negative z-score: %.1f\n', REFINE.max_negative_zscore);
fprintf('  Max rise/decay ratio: %.1f\n', REFINE.max_rise_decay_ratio);
fprintf('  Max single-frame jump: %.2f dF/F\n', REFINE.max_single_frame_jump);
fprintf('  Max plateau duration: %.1f s\n', REFINE.max_plateau_duration_s);
fprintf('  Decay tau range: %.2f - %.1f s\n', REFINE.min_decay_tau_s, REFINE.max_decay_tau_s);
fprintf('\n');

fprintf('═══ INHERITED CRITERIA ═══\n');
fprintf('  Max pairwise correlation: %.2f\n', REFINE.max_pairwise_corr);
fprintf('  Min center distance: %d px\n', REFINE.min_center_distance_px);
fprintf('  Min event count: %d\n', REFINE.min_event_count);
fprintf('  Max SNR: %.1f\n', REFINE.max_snr);
fprintf('  Require baseline return: %s\n', mat2str(REFINE.require_baseline_return));
fprintf('\n');

% Create summary output directory
if ~exist(SUMMARY_OUTPUT_DIR, 'dir')
    mkdir(SUMMARY_OUTPUT_DIR);
    fprintf('Created summary directory: %s\n', SUMMARY_OUTPUT_DIR);
end

% Master summary table
all_summary = table();

%% ════════════════════════════════════════════════════════════════════════
%  PROCESS EACH CONDITION
% ════════════════════════════════════════════════════════════════════════

for c = 1:length(CONDITIONS)
    cond = CONDITIONS(c);
    
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════════════\n');
    fprintf('  CONDITION: %s\n', cond.name);
    fprintf('════════════════════════════════════════════════════════════════════\n');
    
    % Check if input directory exists
    if ~exist(cond.roi_analysis_dir, 'dir')
        fprintf('  [SKIP] Directory not found: %s\n', cond.roi_analysis_dir);
        continue;
    end
    
    % Create condition output folder with subfolders
    cond_output = cond.output_dir;
    if ~exist(cond_output, 'dir')
        mkdir(cond_output);
    end
    fprintf('  Output: %s\n', cond_output);
    
    % Create organized subfolders
    if REFINE.organize_subfolders
        data_dir = fullfile(cond_output, 'data');
        figures_dir = fullfile(cond_output, 'figures');
        logs_dir = fullfile(cond_output, 'logs');
        rejected_dir = fullfile(cond_output, 'rejected_examples');
        
        if ~exist(data_dir, 'dir'), mkdir(data_dir); end
        if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
        if ~exist(logs_dir, 'dir'), mkdir(logs_dir); end
        if ~exist(rejected_dir, 'dir'), mkdir(rejected_dir); end
    else
        data_dir = cond_output;
        figures_dir = cond_output;
        logs_dir = cond_output;
        rejected_dir = cond_output;
    end
    
    % Find ROI .mat files
    mat_files = dir(fullfile(cond.roi_analysis_dir, '*_rois.mat'));
    
    if isempty(mat_files)
        fprintf('  [SKIP] No ROI files found\n');
        continue;
    end
    
    fprintf('  Found %d ROI files\n\n', length(mat_files));
    
    % Process each file
    for f = 1:length(mat_files)
        mat_path = fullfile(mat_files(f).folder, mat_files(f).name);
        [~, fname, ~] = fileparts(mat_files(f).name);
        sample_name = strrep(fname, '_rois', '');
        
        fprintf('─── %s ───\n', sample_name);
        
        % Check if already processed
        if REFINE.skip_existing
            refined_mat_path = fullfile(data_dir, [sample_name '_refined.mat']);
            if exist(refined_mat_path, 'file')
                fprintf('  [SKIP] Already refined (set skip_existing=false to rerun)\n\n');
                continue;
            end
        end
        
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
            fprintf('  [ERROR] Failed to load: %s\n', ME.message);
            continue;
        end
        
        % Store original for visualization
        roi_data_original = roi_data;
        
        % Get ROI count
        n_rois_original = size(roi_data.roi_centers, 1);
        fprintf('  Original ROIs: %d\n', n_rois_original);
        
        if n_rois_original == 0
            fprintf('  [SKIP] No ROIs to refine\n\n');
            continue;
        end
        
        % Initialize removal tracking
        remove_idx = false(n_rois_original, 1);
        removal_reasons = cell(n_rois_original, 1);
        removal_reasons(:) = {''};
        removal_category = cell(n_rois_original, 1);
        removal_category(:) = {''};
        
        % Get traces (ensure correct orientation: [time x ROIs])
        dFF = double(roi_data.dFF);
        if size(dFF, 1) == n_rois_original && size(dFF, 2) ~= n_rois_original
            dFF = dFF';  % Transpose to [time x ROIs]
        end
        
        % Also get smoothed traces if available
        if isfield(roi_data, 'dFF_smooth')
            dFF_smooth = double(roi_data.dFF_smooth);
            if size(dFF_smooth, 1) == n_rois_original && size(dFF_smooth, 2) ~= n_rois_original
                dFF_smooth = dFF_smooth';
            end
        else
            dFF_smooth = dFF;
        end
        
        n_frames = size(dFF, 1);
        
        % Get brain mask
        if isfield(roi_data, 'brain_mask')
            brain_mask = roi_data.brain_mask;
        else
            brain_mask = true(size(roi_data.meanImg));
        end
        
        % Get SNR values
        if isfield(roi_data, 'snr')
            snr_vals = roi_data.snr(:);
        else
            snr_vals = ones(n_rois_original, 1) * 10;
        end
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 1: EXTREME SNR (inherited from v2.3)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [1/10] Checking extreme SNR (>%.1f)...\n', REFINE.max_snr);
        n_removed_snr = 0;
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            if snr_vals(i) > REFINE.max_snr
                remove_idx(i) = true;
                removal_reasons{i} = sprintf('SNR too high (%.1f > %.1f)', snr_vals(i), REFINE.max_snr);
                removal_category{i} = 'extreme_snr';
                n_removed_snr = n_removed_snr + 1;
                continue;
            end
            
            % Secondary check: sustained high signal [v3.2: relaxed 5x→8x, 15%→30%]
            trace = dFF(:, i);
            baseline_region = trace(trace < prctile(trace, 75));
            if isempty(baseline_region), baseline_region = trace; end
            noise_level = 1.4826 * median(abs(baseline_region - median(baseline_region)));
            if noise_level < 0.001, noise_level = 0.001; end
            
            high_signal_fraction = mean(trace > 8 * noise_level);
            if high_signal_fraction > 0.30
                remove_idx(i) = true;
                removal_reasons{i} = sprintf('Sustained high signal (%.0f%% > 8x noise)', high_signal_fraction * 100);
                removal_category{i} = 'extreme_snr';
                n_removed_snr = n_removed_snr + 1;
            end
        end
        fprintf('    Removed: %d\n', n_removed_snr);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 2: EDGE CONTACT (v3.0 NEW)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [2/10] Checking edge contact...\n');
        n_removed_edge = 0;
        
        brain_edge = bwperim(brain_mask);
        
        if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
            for i = 1:n_rois_original
                if remove_idx(i), continue; end
                
                mask = roi_data.roi_masks(:,:,i);
                roi_perim = bwperim(mask);
                
                edge_overlap = roi_perim & brain_edge;
                edge_contact_fraction = sum(edge_overlap(:)) / (sum(roi_perim(:)) + eps);
                
                if edge_contact_fraction > REFINE.max_edge_contact_fraction
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Edge contact (%.0f%% > %.0f%%)', ...
                        edge_contact_fraction * 100, REFINE.max_edge_contact_fraction * 100);
                    removal_category{i} = 'edge_contact';
                    n_removed_edge = n_removed_edge + 1;
                end
            end
        end
        fprintf('    Removed: %d\n', n_removed_edge);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 3: SHAPE METRICS (v3.0 ENHANCED)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [3/10] Checking shape metrics...\n');
        n_removed_shape = 0;
        
        if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
            for i = 1:n_rois_original
                if remove_idx(i), continue; end
                
                mask = roi_data.roi_masks(:,:,i);
                props = regionprops(mask, 'Area', 'Perimeter', 'Solidity', ...
                    'Eccentricity', 'BoundingBox', 'Extent', 'MajorAxisLength', 'MinorAxisLength');
                
                if isempty(props), continue; end
                p = props(1);
                
                % Compactness: 4*pi*area/perimeter^2
                compactness = 4 * pi * p.Area / (p.Perimeter^2 + eps);
                
                % Aspect ratio
                if p.MinorAxisLength > 0
                    aspect_ratio = p.MajorAxisLength / p.MinorAxisLength;
                else
                    aspect_ratio = 999;
                end
                
                % Check all shape criteria
                if compactness < REFINE.min_compactness
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Low compactness (%.2f < %.2f)', compactness, REFINE.min_compactness);
                    removal_category{i} = 'shape_compactness';
                    n_removed_shape = n_removed_shape + 1;
                    continue;
                end
                
                if p.Extent < REFINE.min_extent
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Low extent (%.2f < %.2f) - stringy', p.Extent, REFINE.min_extent);
                    removal_category{i} = 'shape_extent';
                    n_removed_shape = n_removed_shape + 1;
                    continue;
                end
                
                if p.Extent > REFINE.max_extent
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('High extent (%.2f > %.2f) - rectangular', p.Extent, REFINE.max_extent);
                    removal_category{i} = 'shape_extent';
                    n_removed_shape = n_removed_shape + 1;
                    continue;
                end
                
                if aspect_ratio > REFINE.max_aspect_ratio
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('High aspect ratio (%.1f > %.1f)', aspect_ratio, REFINE.max_aspect_ratio);
                    removal_category{i} = 'shape_aspect';
                    n_removed_shape = n_removed_shape + 1;
                    continue;
                end
                
                if p.Solidity < REFINE.min_solidity
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Low solidity (%.2f < %.2f)', p.Solidity, REFINE.min_solidity);
                    removal_category{i} = 'shape_solidity';
                    n_removed_shape = n_removed_shape + 1;
                    continue;
                end
                
                if p.Eccentricity > REFINE.max_eccentricity
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('High eccentricity (%.2f > %.2f)', p.Eccentricity, REFINE.max_eccentricity);
                    removal_category{i} = 'shape_eccentricity';
                    n_removed_shape = n_removed_shape + 1;
                end
            end
        end
        fprintf('    Removed: %d\n', n_removed_shape);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 4: NEGATIVE TRANSIENTS (v3.0 NEW)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [4/10] Checking for negative transients...\n');
        n_removed_negative = 0;
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            trace = dFF(:, i);
            
            baseline_region = trace(trace < prctile(trace, 50));
            if length(baseline_region) < 10, baseline_region = trace; end
            baseline = median(baseline_region);
            noise = 1.4826 * median(abs(baseline_region - baseline));
            if noise < 0.001, noise = 0.001; end
            
            min_zscore = (min(trace) - baseline) / noise;
            
            if min_zscore < REFINE.max_negative_zscore
                remove_idx(i) = true;
                removal_reasons{i} = sprintf('Negative transient (z=%.1f < %.1f)', min_zscore, REFINE.max_negative_zscore);
                removal_category{i} = 'kinetics_negative';
                n_removed_negative = n_removed_negative + 1;
            end
        end
        fprintf('    Removed: %d\n', n_removed_negative);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 5: STEP ARTIFACTS (v3.0 NEW)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [5/10] Checking for step artifacts...\n');
        n_removed_step = 0;
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            trace = dFF(:, i);
            frame_diff = abs(diff(trace));
            large_jumps = find(frame_diff > REFINE.max_single_frame_jump);
            
            n_bad_jumps = 0;
            for j = 1:length(large_jumps)
                jump_idx = large_jumps(j);
                
                if jump_idx < n_frames - round(0.5 * FRAME_RATE)
                    post_jump = trace(jump_idx+1 : min(jump_idx + round(0.5*FRAME_RATE), n_frames));
                    pre_jump = trace(max(1, jump_idx-5) : jump_idx);
                    
                    if mean(post_jump) > mean(pre_jump) + 0.8 * frame_diff(jump_idx)
                        n_bad_jumps = n_bad_jumps + 1;
                    end
                end
            end
            
            if n_bad_jumps >= 2
                remove_idx(i) = true;
                removal_reasons{i} = sprintf('Step artifacts (%d jumps stay elevated)', n_bad_jumps);
                removal_category{i} = 'kinetics_step';
                n_removed_step = n_removed_step + 1;
            end
        end
        fprintf('    Removed: %d\n', n_removed_step);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 6: RISE/DECAY KINETICS (v3.0 NEW)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [6/10] Checking rise/decay kinetics...\n');
        n_removed_kinetics = 0;
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            trace = dFF_smooth(:, i);
            
            try
                [pks, locs] = findpeaks(trace, 'MinPeakProminence', 0.05, ...
                    'MinPeakDistance', round(FRAME_RATE * 0.5));
            catch
                pks = []; locs = [];
            end
            
            if length(pks) < 2, continue; end
            
            [~, pk_order] = sort(pks, 'descend');
            n_analyze = min(5, length(pks));
            
            rise_times = [];
            decay_times = [];
            
            for p = 1:n_analyze
                pk_idx = locs(pk_order(p));
                pk_val = pks(pk_order(p));
                
                search_start = max(1, pk_idx - round(2 * FRAME_RATE));
                pre_region = trace(search_start:pk_idx);
                [pre_min, pre_min_idx] = min(pre_region);
                rise_start = search_start + pre_min_idx - 1;
                
                search_end = min(n_frames, pk_idx + round(3 * FRAME_RATE));
                post_region = trace(pk_idx:search_end);
                decay_target = pre_min + 0.37 * (pk_val - pre_min);
                decay_idx = find(post_region < decay_target, 1, 'first');
                
                if ~isempty(decay_idx) && decay_idx > 1
                    rise_time = (pk_idx - rise_start) / FRAME_RATE;
                    decay_time = decay_idx / FRAME_RATE;
                    
                    if rise_time > 0.01 && decay_time > 0.01
                        rise_times(end+1) = rise_time;
                        decay_times(end+1) = decay_time;
                    end
                end
            end
            
            if ~isempty(rise_times) && ~isempty(decay_times)
                median_rise = median(rise_times);
                median_decay = median(decay_times);
                rise_decay_ratio = median_rise / (median_decay + eps);
                
                if rise_decay_ratio > REFINE.max_rise_decay_ratio
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Abnormal kinetics (rise/decay=%.2f)', rise_decay_ratio);
                    removal_category{i} = 'kinetics_ratio';
                    n_removed_kinetics = n_removed_kinetics + 1;
                    continue;
                end
                
                if median_decay < REFINE.min_decay_tau_s
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Decay too fast (tau=%.3fs)', median_decay);
                    removal_category{i} = 'kinetics_decay';
                    n_removed_kinetics = n_removed_kinetics + 1;
                elseif median_decay > REFINE.max_decay_tau_s
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('Decay too slow (tau=%.1fs)', median_decay);
                    removal_category{i} = 'kinetics_decay';
                    n_removed_kinetics = n_removed_kinetics + 1;
                end
            end
        end
        fprintf('    Removed: %d\n', n_removed_kinetics);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 7: SUSTAINED PLATEAUS (v3.0 NEW)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [7/10] Checking for sustained plateaus...\n');
        n_removed_plateau = 0;
        
        plateau_frames = round(REFINE.max_plateau_duration_s * FRAME_RATE);
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            trace = dFF(:, i);
            baseline = prctile(trace, 10);
            peak_val = prctile(trace, 99);
            
            if peak_val - baseline < 0.05, continue; end
            
            elevated_thresh = baseline + REFINE.plateau_threshold * (peak_val - baseline);
            elevated = trace > elevated_thresh;
            
            d = diff([0; elevated; 0]);
            starts = find(d == 1);
            ends = find(d == -1) - 1;
            run_lengths = ends - starts + 1;
            
            if ~isempty(run_lengths) && max(run_lengths) > plateau_frames
                remove_idx(i) = true;
                removal_reasons{i} = sprintf('Sustained plateau (%.1fs > %.1fs)', ...
                    max(run_lengths) / FRAME_RATE, REFINE.max_plateau_duration_s);
                removal_category{i} = 'kinetics_plateau';
                n_removed_plateau = n_removed_plateau + 1;
            end
        end
        fprintf('    Removed: %d\n', n_removed_plateau);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 8: CORRELATION (inherited from v2.3)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [8/10] Checking pairwise correlations...\n');
        n_removed_corr = 0;
        
        corr_matrix = corrcoef(dFF);
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            for j = (i+1):n_rois_original
                if remove_idx(j), continue; end
                
                if abs(corr_matrix(i, j)) > REFINE.max_pairwise_corr
                    if REFINE.corr_keep_higher_snr
                        if snr_vals(i) >= snr_vals(j)
                            remove_idx(j) = true;
                            removal_reasons{j} = sprintf('Correlated with ROI %d (r=%.2f)', i, corr_matrix(i,j));
                            removal_category{j} = 'correlation';
                        else
                            remove_idx(i) = true;
                            removal_reasons{i} = sprintf('Correlated with ROI %d (r=%.2f)', j, corr_matrix(i,j));
                            removal_category{i} = 'correlation';
                        end
                    else
                        remove_idx(j) = true;
                        removal_reasons{j} = sprintf('Correlated with ROI %d (r=%.2f)', i, corr_matrix(i,j));
                        removal_category{j} = 'correlation';
                    end
                    n_removed_corr = n_removed_corr + 1;
                end
            end
        end
        fprintf('    Removed: %d\n', n_removed_corr);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 9: SPATIAL PROXIMITY (inherited from v2.3)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [9/10] Checking spatial proximity...\n');
        n_removed_spatial = 0;
        
        centers = roi_data.roi_centers;
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            for j = (i+1):n_rois_original
                if remove_idx(j), continue; end
                
                dist = sqrt((centers(i,1) - centers(j,1))^2 + (centers(i,2) - centers(j,2))^2);
                
                if dist < REFINE.min_center_distance_px
                    if snr_vals(i) >= snr_vals(j)
                        remove_idx(j) = true;
                        removal_reasons{j} = sprintf('Too close to ROI %d (%.1f px)', i, dist);
                        removal_category{j} = 'spatial';
                    else
                        remove_idx(i) = true;
                        removal_reasons{i} = sprintf('Too close to ROI %d (%.1f px)', j, dist);
                        removal_category{i} = 'spatial';
                    end
                    n_removed_spatial = n_removed_spatial + 1;
                end
            end
        end
        fprintf('    Removed: %d\n', n_removed_spatial);
        
        % ═══════════════════════════════════════════════════════════════════
        % CRITERION 10: BASELINE DRIFT & ACTIVITY (inherited from v2.3)
        % ═══════════════════════════════════════════════════════════════════
        fprintf('  [10/10] Checking activity & baseline...\n');
        n_removed_activity = 0;
        n_removed_baseline = 0;
        
        baseline_window_frames = round(REFINE.baseline_window_s * FRAME_RATE);
        baseline_window_frames = min(baseline_window_frames, floor(n_frames / 4));
        
        for i = 1:n_rois_original
            if remove_idx(i), continue; end
            
            % Activity check
            if isfield(roi_data, 'event_count')
                if roi_data.event_count(i) < REFINE.min_event_count
                    remove_idx(i) = true;
                    removal_reasons{i} = sprintf('No events (%d < %d)', roi_data.event_count(i), REFINE.min_event_count);
                    removal_category{i} = 'inactive';
                    n_removed_activity = n_removed_activity + 1;
                    continue;
                end
            end
            
            % Baseline drift check
            if REFINE.require_baseline_return
                trace = dFF(:, i);
                
                start_window = trace(1:baseline_window_frames);
                end_idx_start = max(1, n_frames - baseline_window_frames + 1);
                end_window = trace(end_idx_start:n_frames);
                
                start_trimmed = start_window(start_window <= prctile(start_window, 90));
                end_trimmed = end_window(end_window <= prctile(end_window, 90));
                
                if isempty(start_trimmed), start_trimmed = start_window; end
                if isempty(end_trimmed), end_trimmed = end_window; end
                
                start_baseline = median(start_trimmed);
                end_baseline = median(end_trimmed);
                drift = end_baseline - start_baseline;
                
                if abs(drift) > REFINE.baseline_return_threshold
                    elevated_threshold = start_baseline + REFINE.baseline_return_threshold * 0.5;
                    fraction_elevated = mean(end_window > elevated_threshold);
                    
                    if fraction_elevated > 0.5
                        remove_idx(i) = true;
                        removal_reasons{i} = sprintf('Baseline drift (%.2f dF/F, %.0f%% elevated)', drift, fraction_elevated * 100);
                        removal_category{i} = 'baseline_drift';
                        n_removed_baseline = n_removed_baseline + 1;
                    end
                end
            end
        end
        fprintf('    Removed (inactive): %d, (baseline): %d\n', n_removed_activity, n_removed_baseline);
        
        % ═══════════════════════════════════════════════════════════════════
        % CREATE REFINED DATA
        % ═══════════════════════════════════════════════════════════════════
        keep_idx = ~remove_idx;
        n_rois_refined = sum(keep_idx);
        
        fprintf('  ────────────────────────────────────────\n');
        fprintf('  REFINEMENT SUMMARY:\n');
        fprintf('    Original: %d ROIs\n', n_rois_original);
        fprintf('    Removed:  %d ROIs (%.1f%%)\n', sum(remove_idx), 100*sum(remove_idx)/n_rois_original);
        fprintf('      - Extreme SNR: %d\n', n_removed_snr);
        fprintf('      - Edge contact: %d\n', n_removed_edge);
        fprintf('      - Shape: %d\n', n_removed_shape);
        fprintf('      - Negative transients: %d\n', n_removed_negative);
        fprintf('      - Step artifacts: %d\n', n_removed_step);
        fprintf('      - Kinetics: %d\n', n_removed_kinetics);
        fprintf('      - Plateau: %d\n', n_removed_plateau);
        fprintf('      - Correlated: %d\n', n_removed_corr);
        fprintf('      - Spatial: %d\n', n_removed_spatial);
        fprintf('      - Inactive: %d\n', n_removed_activity);
        fprintf('      - Baseline: %d\n', n_removed_baseline);
        fprintf('    Refined:  %d ROIs\n', n_rois_refined);
        fprintf('  ────────────────────────────────────────\n');
        
        % Build refined roi_data structure
        refined_data = struct();
        
        % Copy and filter arrays
        refined_data.roi_centers = roi_data.roi_centers(keep_idx, :);
        
        if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
            refined_data.roi_masks = roi_data.roi_masks(:, :, keep_idx);
        end
        
        % Handle dFF orientation
        original_dFF = roi_data.dFF;
        if size(original_dFF, 1) == n_rois_original
            refined_data.dFF = original_dFF(keep_idx, :);
        else
            refined_data.dFF = original_dFF(:, keep_idx);
        end
        
        if isfield(roi_data, 'dFF_smooth')
            original_smooth = roi_data.dFF_smooth;
            if size(original_smooth, 1) == n_rois_original
                refined_data.dFF_smooth = original_smooth(keep_idx, :);
            else
                refined_data.dFF_smooth = original_smooth(:, keep_idx);
            end
        end
        
        if isfield(roi_data, 'F')
            original_F = roi_data.F;
            if size(original_F, 1) == n_rois_original
                refined_data.F = original_F(keep_idx, :);
            else
                refined_data.F = original_F(:, keep_idx);
            end
        end
        
        if isfield(roi_data, 'Fneu')
            original_Fneu = roi_data.Fneu;
            if size(original_Fneu, 1) == n_rois_original
                refined_data.Fneu = original_Fneu(keep_idx, :);
            else
                refined_data.Fneu = original_Fneu(:, keep_idx);
            end
        end
        
        % v3.2/v4.5 compatibility: Handle dFF_raw (pre-regression traces)
        if isfield(roi_data, 'dFF_raw')
            original_dFF_raw = roi_data.dFF_raw;
            if size(original_dFF_raw, 1) == n_rois_original
                refined_data.dFF_raw = original_dFF_raw(keep_idx, :);
            else
                refined_data.dFF_raw = original_dFF_raw(:, keep_idx);
            end
        end
        
        % Filter scalar metrics
        if isfield(roi_data, 'snr')
            refined_data.snr = roi_data.snr(keep_idx);
        end
        if isfield(roi_data, 'event_rate')
            refined_data.event_rate = roi_data.event_rate(keep_idx);
        end
        if isfield(roi_data, 'event_count')
            refined_data.event_count = roi_data.event_count(keep_idx);
        end
        if isfield(roi_data, 'event_amplitude')
            refined_data.event_amplitude = roi_data.event_amplitude(keep_idx);
        end
        
        % Copy non-ROI-specific fields
        if isfield(roi_data, 'meanImg')
            refined_data.meanImg = roi_data.meanImg;
        end
        if isfield(roi_data, 'Cn')
            refined_data.Cn = roi_data.Cn;
        end
        if isfield(roi_data, 'PNR')
            refined_data.PNR = roi_data.PNR;
        end
        if isfield(roi_data, 'brain_mask')
            refined_data.brain_mask = roi_data.brain_mask;
        end
        if isfield(roi_data, 'detection_info')
            refined_data.detection_info = roi_data.detection_info;
        end
        
        % v3.2/v4.5 compatibility: Copy quality_metrics (global correlation, neuropil severity, etc.)
        if isfield(roi_data, 'quality_metrics')
            refined_data.quality_metrics = roi_data.quality_metrics;
            % Filter per-ROI quality metrics if present
            if isfield(roi_data.quality_metrics, 'global_corr_per_roi_before')
                refined_data.quality_metrics.global_corr_per_roi_before = ...
                    roi_data.quality_metrics.global_corr_per_roi_before(keep_idx);
            end
            if isfield(roi_data.quality_metrics, 'global_corr_per_roi_after')
                refined_data.quality_metrics.global_corr_per_roi_after = ...
                    roi_data.quality_metrics.global_corr_per_roi_after(keep_idx);
            end
        end
        
        % Add original index mapping
        original_indices = find(keep_idx);
        refined_data.original_roi_idx = original_indices;
        
        % Recompute synchrony for refined ROIs
        if n_rois_refined >= 2
            dFF_refined = dFF(:, keep_idx);
            corr_refined = corrcoef(dFF_refined);
            upper_tri = triu(true(n_rois_refined), 1);
            refined_data.synchrony_score = mean(abs(corr_refined(upper_tri)), 'omitnan');
            refined_data.synchrony_matrix = corr_refined;
        else
            refined_data.synchrony_score = NaN;
            refined_data.synchrony_matrix = [];
        end
        
        % Add provenance
        refined_data.refinement = struct();
        refined_data.refinement.date = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
        refined_data.refinement.version = '3.2';
        refined_data.refinement.n_original = n_rois_original;
        refined_data.refinement.n_refined = n_rois_refined;
        refined_data.refinement.n_removed = sum(remove_idx);
        refined_data.refinement.criteria = REFINE;
        refined_data.refinement.removal_counts = struct(...
            'extreme_snr', n_removed_snr, ...
            'edge_contact', n_removed_edge, ...
            'shape', n_removed_shape, ...
            'negative_transients', n_removed_negative, ...
            'step_artifacts', n_removed_step, ...
            'kinetics', n_removed_kinetics, ...
            'plateau', n_removed_plateau, ...
            'correlation', n_removed_corr, ...
            'spatial', n_removed_spatial, ...
            'inactive', n_removed_activity, ...
            'baseline_drift', n_removed_baseline);
        
        % Copy original provenance
        if isfield(roi_data, 'provenance')
            refined_data.original_provenance = roi_data.provenance;
            if isfield(roi_data.provenance, 'motion_quality')
                motion_quality = char(roi_data.provenance.motion_quality);
            else
                motion_quality = 'unknown';
            end
        else
            motion_quality = 'unknown';
        end
        
        % ═══════════════════════════════════════════════════════════════════
        % SAVE REFINED DATA
        % ═══════════════════════════════════════════════════════════════════
        
        refined_mat_path = fullfile(data_dir, [sample_name '_refined.mat']);
        roi_data = refined_data;
        save(refined_mat_path, 'roi_data', '-v7.3');
        fprintf('  Saved: %s\n', refined_mat_path);
        
        % ═══════════════════════════════════════════════════════════════════
        % SAVE REFINEMENT LOG
        % ═══════════════════════════════════════════════════════════════════
        
        log_path = fullfile(logs_dir, [sample_name '_refinement_log.csv']);
        
        log_table = table();
        log_table.ROI_Original = (1:n_rois_original)';
        log_table.Kept = keep_idx;
        log_table.Category = removal_category;
        log_table.Removal_Reason = removal_reasons;
        log_table.SNR = snr_vals;
        
        writetable(log_table, log_path);
        fprintf('  Saved: %s\n', log_path);
        
        % ═══════════════════════════════════════════════════════════════════
        % GENERATE VISUALIZATIONS
        % ═══════════════════════════════════════════════════════════════════
        
        if REFINE.generate_gallery && n_rois_refined > 0
            fprintf('  Generating visualizations...\n');
            
            % Generate rejection examples figure (v3.0 NEW)
            if REFINE.generate_rejection_examples && sum(remove_idx) > 0
                generate_rejection_examples_v3(roi_data_original, remove_idx, removal_reasons, ...
                    removal_category, dFF, rejected_dir, sample_name, FRAME_RATE);
            end
            
            % Generate summary figure showing removed vs kept ROIs
            generate_refinement_summary_figure_v3(roi_data_original, keep_idx, removal_reasons, ...
                removal_category, figures_dir, sample_name);
            
            % Generate improved ROI gallery
            generate_roi_gallery_v3(refined_data, dFF(:, keep_idx), ...
                figures_dir, sample_name, FRAME_RATE, REFINE);
            
            % Generate analysis figures
            if REFINE.generate_analysis_figures
                generate_refined_analysis_figure_v3(refined_data, figures_dir, sample_name, FRAME_RATE);
            end
        end
        
        % ═══════════════════════════════════════════════════════════════════
        % ADD TO SUMMARY TABLE
        % ═══════════════════════════════════════════════════════════════════
        
        row = table();
        row.filename = {sample_name};
        row.condition = {cond.name};
        row.motion_quality = {motion_quality};
        row.n_original = n_rois_original;
        row.n_refined = n_rois_refined;
        row.n_removed = sum(remove_idx);
        row.pct_removed = 100 * sum(remove_idx) / n_rois_original;
        row.removed_snr = n_removed_snr;
        row.removed_edge = n_removed_edge;
        row.removed_shape = n_removed_shape;
        row.removed_negative = n_removed_negative;
        row.removed_step = n_removed_step;
        row.removed_kinetics = n_removed_kinetics;
        row.removed_plateau = n_removed_plateau;
        row.removed_corr = n_removed_corr;
        row.removed_spatial = n_removed_spatial;
        row.removed_inactive = n_removed_activity;
        row.removed_baseline = n_removed_baseline;
        
        if n_rois_refined > 0
            row.mean_snr = mean(refined_data.snr);
            row.mean_event_rate = mean(refined_data.event_rate);
            row.synchrony_score = refined_data.synchrony_score;
        else
            row.mean_snr = NaN;
            row.mean_event_rate = NaN;
            row.synchrony_score = NaN;
        end
        
        all_summary = [all_summary; row];
        
        fprintf('\n');
    end
end

%% ════════════════════════════════════════════════════════════════════════
%  SAVE MASTER SUMMARY
% ════════════════════════════════════════════════════════════════════════

if ~isempty(all_summary)
    % Save summary CSV
    summary_path = fullfile(SUMMARY_OUTPUT_DIR, 'refinement_summary.csv');
    writetable(all_summary, summary_path);
    fprintf('\n✓ Saved summary: %s\n', summary_path);
    
    % Create group comparison
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════════════\n');
    fprintf('  GROUP COMPARISON SUMMARY\n');
    fprintf('════════════════════════════════════════════════════════════════════\n');
    
    conditions_in_data = unique(all_summary.condition);
    group_comparison = table();
    
    for i = 1:length(conditions_in_data)
        cond_name = conditions_in_data{i};
        cond_rows = strcmp(all_summary.condition, cond_name);
        
        row = table();
        row.condition = {cond_name};
        row.n_files = sum(cond_rows);
        row.total_rois_original = sum(all_summary.n_original(cond_rows));
        row.total_rois_refined = sum(all_summary.n_refined(cond_rows));
        row.mean_rois_per_file = mean(all_summary.n_refined(cond_rows));
        row.std_rois_per_file = std(all_summary.n_refined(cond_rows));
        row.mean_pct_removed = mean(all_summary.pct_removed(cond_rows));
        row.total_edge_removed = sum(all_summary.removed_edge(cond_rows));
        row.total_shape_removed = sum(all_summary.removed_shape(cond_rows));
        row.total_kinetics_removed = sum(all_summary.removed_kinetics(cond_rows) + ...
            all_summary.removed_negative(cond_rows) + all_summary.removed_step(cond_rows) + ...
            all_summary.removed_plateau(cond_rows));
        row.mean_snr = mean(all_summary.mean_snr(cond_rows), 'omitnan');
        row.std_snr = std(all_summary.mean_snr(cond_rows), 'omitnan');
        row.mean_event_rate = mean(all_summary.mean_event_rate(cond_rows), 'omitnan');
        row.std_event_rate = std(all_summary.mean_event_rate(cond_rows), 'omitnan');
        row.mean_synchrony = mean(all_summary.synchrony_score(cond_rows), 'omitnan');
        
        group_comparison = [group_comparison; row];
        
        fprintf('  %s:\n', cond_name);
        fprintf('    Files: %d\n', row.n_files);
        fprintf('    Total ROIs: %d → %d (%.1f%% removed)\n', ...
            row.total_rois_original, row.total_rois_refined, row.mean_pct_removed);
        fprintf('    v3.0 removals - Edge: %d, Shape: %d, Kinetics: %d\n', ...
            row.total_edge_removed, row.total_shape_removed, row.total_kinetics_removed);
        fprintf('    Mean ROIs/file: %.1f ± %.1f\n', row.mean_rois_per_file, row.std_rois_per_file);
        fprintf('    Mean SNR: %.2f ± %.2f\n', row.mean_snr, row.std_snr);
        fprintf('    Mean Event Rate: %.3f ± %.3f Hz\n', row.mean_event_rate, row.std_event_rate);
        fprintf('    Mean Synchrony: %.2f\n', row.mean_synchrony);
        fprintf('\n');
    end
    
    group_path = fullfile(SUMMARY_OUTPUT_DIR, 'group_comparison.csv');
    writetable(group_comparison, group_path);
    fprintf('✓ Saved group comparison: %s\n', group_path);
    
    % Generate README
    generate_readme_v3(SUMMARY_OUTPUT_DIR, all_summary, REFINE, CONDITIONS);
    fprintf('✓ Saved README.txt\n');
end

%% ════════════════════════════════════════════════════════════════════════
%  COMPLETION
% ════════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 ROI REFINEMENT v3.2 COMPLETE                     ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Output locations:\n');
for c = 1:length(CONDITIONS)
    fprintf('  %s: %s\n', CONDITIONS(c).name, CONDITIONS(c).output_dir);
end
fprintf('  Summaries: %s\n', SUMMARY_OUTPUT_DIR);
fprintf('\n');
fprintf('v3.0 NEW features:\n');
fprintf('  - Edge contact removal (ROIs touching brain boundary)\n');
fprintf('  - Shape validation (compactness, extent, aspect ratio)\n');
fprintf('  - Kinetics validation (rise/decay ratio, plateaus)\n');
fprintf('  - Negative transient detection\n');
fprintf('  - Step artifact detection\n');
fprintf('  - Rejection examples figure (see rejected_examples/ folder)\n');
fprintf('\n');

end  % end of main function ROI_REFINEMENT_v3_2

%% ════════════════════════════════════════════════════════════════════════
%  CONFIG HELPERS
% ════════════════════════════════════════════════════════════════════════

function cfg = merge_with_defaults(user_cfg)
% Merge a user-provided config struct with defaults.
cfg = roi_refinement_config();
default_fields = fieldnames(cfg);
user_fields = fieldnames(user_cfg);

for i = 1:numel(user_fields)
    if ismember(user_fields{i}, default_fields)
        cfg.(user_fields{i}) = user_cfg.(user_fields{i});
    else
        warning('roi_refinement:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See roi_refinement_config() for valid fields.', ...
            user_fields{i});
    end
end
end

function conditions = prompt_for_conditions()
% Prompt the user to define conditions via folder pickers.
if ~usejava('desktop')
    error(['No conditions defined and no GUI available.\n' ...
           'Usage:\n' ...
           '  cfg = roi_refinement_config();\n' ...
           '  cfg.conditions(1).name = ''Control'';\n' ...
           '  cfg.conditions(1).roi_analysis_dir = ''/path/to/roi_analysis'';\n' ...
           '  cfg.conditions(1).output_dir = ''/path/to/refined'';\n' ...
           '  ROI_REFINEMENT_v3_2(cfg);\n\n' ...
           'See also: roi_refinement_config']);
end

conditions = struct('name', {}, 'roi_analysis_dir', {}, 'output_dir', {});
condIdx = 0;

while true
    condIdx = condIdx + 1;
    answer = questdlg(sprintf('Add condition %d?', condIdx), ...
        'Add Condition', 'Yes', 'No (done)', 'Yes');
    if ~strcmp(answer, 'Yes')
        break;
    end
    
    name = inputdlg(sprintf('Condition %d name:', condIdx), 'Condition Name', 1, {sprintf('Condition_%d', condIdx)});
    if isempty(name), break; end
    
    roiDir = uigetdir('', sprintf('Select ROI analysis folder for %s', name{1}));
    if roiDir == 0, break; end
    
    outDir = uigetdir('', sprintf('Select output folder for %s', name{1}));
    if outDir == 0, break; end
    
    conditions(condIdx).name = name{1};
    conditions(condIdx).roi_analysis_dir = roiDir;
    conditions(condIdx).output_dir = outDir;
end

if isempty(conditions)
    error('No conditions defined. Exiting.');
end
end

%% ════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
% ════════════════════════════════════════════════════════════════════════

function generate_rejection_examples_v3(roi_data, remove_idx, reasons, categories, dFF, output_dir, sample_name, frame_rate)
% v3.0 NEW: Generate figure showing examples of WHY ROIs were rejected
% Groups by rejection category

removed_cats = categories(remove_idx);
unique_cats = unique(removed_cats(~cellfun(@isempty, removed_cats)));

if isempty(unique_cats)
    return;
end

n_frames = size(dFF, 1);
time_vec = (0:n_frames-1) / frame_rate;

% Create figure showing examples from each category
fig = figure('Position', [50 50 1800 1200], 'Visible', 'off');

n_cats = min(6, length(unique_cats));
n_examples_per_cat = 2;

for cat_idx = 1:n_cats
    cat_name = unique_cats{cat_idx};
    cat_rois = find(remove_idx & strcmp(categories, cat_name));
    
    if isempty(cat_rois), continue; end
    
    n_show = min(n_examples_per_cat, length(cat_rois));
    
    for ex = 1:n_show
        roi_num = cat_rois(ex);
        
        % Show ROI mask
        subplot(n_cats, n_examples_per_cat * 2, (cat_idx-1) * n_examples_per_cat * 2 + (ex-1)*2 + 1);
        
        if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
            mask = roi_data.roi_masks(:,:,roi_num);
            center = roi_data.roi_centers(roi_num, :);
            
            crop_size = 40;
            y_range = max(1, round(center(2))-crop_size):min(size(mask,1), round(center(2))+crop_size);
            x_range = max(1, round(center(1))-crop_size):min(size(mask,2), round(center(1))+crop_size);
            
            if isfield(roi_data, 'meanImg')
                crop_img = roi_data.meanImg(y_range, x_range);
            else
                crop_img = zeros(length(y_range), length(x_range));
            end
            crop_mask = mask(y_range, x_range);
            
            imagesc(crop_img);
            colormap(gca, gray);
            hold on;
            [B, ~] = bwboundaries(crop_mask);
            if ~isempty(B)
                plot(B{1}(:,2), B{1}(:,1), 'r', 'LineWidth', 2);
            end
            axis image off;
        end
        
        title(sprintf('%s\nROI %d', strrep(cat_name, '_', ' '), roi_num), ...
            'FontSize', 8, 'Color', 'r');
        
        % Show trace
        subplot(n_cats, n_examples_per_cat * 2, (cat_idx-1) * n_examples_per_cat * 2 + (ex-1)*2 + 2);
        
        trace = dFF(:, roi_num);
        plot(time_vec, trace, 'r', 'LineWidth', 0.5);
        hold on;
        yline(0, 'k--', 'Alpha', 0.3);
        
        xlabel('Time (s)', 'FontSize', 7);
        ylabel('dF/F', 'FontSize', 7);
        
        reason_str = reasons{roi_num};
        if length(reason_str) > 50
            reason_str = [reason_str(1:47) '...'];
        end
        title(reason_str, 'FontSize', 7, 'Interpreter', 'none');
        
        xlim([0, min(60, time_vec(end))]);
        set(gca, 'FontSize', 7);
    end
end

sgtitle(sprintf('%s - Rejection Examples (v3.0)', sample_name), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

saveas(fig, fullfile(output_dir, [sample_name '_rejection_examples.png']));
close(fig);

fprintf('    ✓ Generated rejection examples figure\n');

end


function generate_refinement_summary_figure_v3(roi_data, keep_idx, removal_reasons, removal_category, output_dir, sample_name)
% Summary figure: kept vs removed ROIs with v3.0 categories

n_rois = length(keep_idx);
n_kept = sum(keep_idx);
n_removed = sum(~keep_idx);

fig = figure('Position', [50 50 1400 1000], 'Visible', 'off');

% Panel 1: ROI map - kept (green) vs removed (red)
subplot(2, 2, 1);
if isfield(roi_data, 'meanImg')
    imagesc(roi_data.meanImg);
    colormap(gca, gray);
else
    imagesc(zeros(100, 100));
end
axis image; hold on;

% Removed ROIs in red
removed_idx = find(~keep_idx);
for i = 1:length(removed_idx)
    r = removed_idx(i);
    if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
        [B, ~] = bwboundaries(roi_data.roi_masks(:,:,r));
        if ~isempty(B)
            plot(B{1}(:,2), B{1}(:,1), 'r', 'LineWidth', 1.5);
        end
    end
    plot(roi_data.roi_centers(r,1), roi_data.roi_centers(r,2), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
end

% Kept ROIs in green
kept_idx_list = find(keep_idx);
for i = 1:length(kept_idx_list)
    r = kept_idx_list(i);
    if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
        [B, ~] = bwboundaries(roi_data.roi_masks(:,:,r));
        if ~isempty(B)
            plot(B{1}(:,2), B{1}(:,1), 'g', 'LineWidth', 1.5);
        end
    end
    plot(roi_data.roi_centers(r,1), roi_data.roi_centers(r,2), 'g+', 'MarkerSize', 6, 'LineWidth', 1.5);
end

title(sprintf('ROI Map: %d kept (green), %d removed (red)', n_kept, n_removed), 'FontSize', 10);

% Panel 2: Pie chart by category
subplot(2, 2, 2);

% Count by category
cat_counts = struct();
cat_counts.extreme_snr = sum(strcmp(removal_category, 'extreme_snr'));
cat_counts.edge_contact = sum(strcmp(removal_category, 'edge_contact'));
cat_counts.shape = sum(contains(removal_category, 'shape'));
cat_counts.kinetics = sum(contains(removal_category, 'kinetics'));
cat_counts.correlation = sum(strcmp(removal_category, 'correlation'));
cat_counts.spatial = sum(strcmp(removal_category, 'spatial'));
cat_counts.inactive = sum(strcmp(removal_category, 'inactive'));
cat_counts.baseline = sum(strcmp(removal_category, 'baseline_drift'));

labels = {};
counts = [];
colors_pie = {};

if cat_counts.edge_contact > 0
    labels{end+1} = sprintf('Edge (%d)', cat_counts.edge_contact);
    counts(end+1) = cat_counts.edge_contact;
    colors_pie{end+1} = [0.8 0.4 0.0];  % Orange
end
if cat_counts.shape > 0
    labels{end+1} = sprintf('Shape (%d)', cat_counts.shape);
    counts(end+1) = cat_counts.shape;
    colors_pie{end+1} = [0.6 0.2 0.8];  % Purple
end
if cat_counts.kinetics > 0
    labels{end+1} = sprintf('Kinetics (%d)', cat_counts.kinetics);
    counts(end+1) = cat_counts.kinetics;
    colors_pie{end+1} = [0.2 0.6 0.8];  % Cyan
end
if cat_counts.extreme_snr > 0
    labels{end+1} = sprintf('SNR (%d)', cat_counts.extreme_snr);
    counts(end+1) = cat_counts.extreme_snr;
    colors_pie{end+1} = [0.9 0.1 0.1];  % Red
end
if cat_counts.correlation > 0
    labels{end+1} = sprintf('Correlated (%d)', cat_counts.correlation);
    counts(end+1) = cat_counts.correlation;
    colors_pie{end+1} = [0.3 0.7 0.3];  % Green
end
if cat_counts.spatial > 0
    labels{end+1} = sprintf('Spatial (%d)', cat_counts.spatial);
    counts(end+1) = cat_counts.spatial;
    colors_pie{end+1} = [0.5 0.5 0.5];  % Gray
end
if cat_counts.inactive > 0
    labels{end+1} = sprintf('Inactive (%d)', cat_counts.inactive);
    counts(end+1) = cat_counts.inactive;
    colors_pie{end+1} = [0.7 0.7 0.2];  % Yellow
end
if cat_counts.baseline > 0
    labels{end+1} = sprintf('Baseline (%d)', cat_counts.baseline);
    counts(end+1) = cat_counts.baseline;
    colors_pie{end+1} = [0.4 0.4 0.8];  % Blue
end

if ~isempty(counts)
    h = pie(counts, labels);
    % Apply colors
    for i = 1:length(counts)
        h(2*i-1).FaceColor = colors_pie{i};
    end
    title('Removal Reasons (v3.0)', 'FontSize', 10);
else
    text(0.5, 0.5, 'No ROIs removed', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12);
    title('Removal Reasons', 'FontSize', 10);
end

% Panel 3: SNR distribution
subplot(2, 2, 3);
if isfield(roi_data, 'snr')
    snr_all = roi_data.snr(:);
    snr_kept = snr_all(keep_idx);
    snr_removed = snr_all(~keep_idx);
    
    hold on;
    if ~isempty(snr_removed)
        histogram(snr_removed, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Removed');
    end
    histogram(snr_kept, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'DisplayName', 'Kept');
    
    xlabel('SNR');
    ylabel('Count');
    title(sprintf('SNR Distribution (kept: %.1f ± %.1f)', mean(snr_kept), std(snr_kept)), 'FontSize', 10);
    legend('Location', 'best');
else
    text(0.5, 0.5, 'No SNR data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
end

% Panel 4: Summary stats
subplot(2, 2, 4);
axis off;

summary_text = {
    sprintf('REFINEMENT SUMMARY v3.0'), ...
    sprintf('%s', sample_name), ...
    '', ...
    sprintf('Original ROIs: %d', n_rois), ...
    sprintf('Kept ROIs: %d (%.1f%%)', n_kept, 100*n_kept/n_rois), ...
    sprintf('Removed ROIs: %d (%.1f%%)', n_removed, 100*n_removed/n_rois), ...
    '', ...
    '─── v3.0 NEW CRITERIA ───', ...
    sprintf('  Edge contact: %d', cat_counts.edge_contact), ...
    sprintf('  Shape (compact/extent/aspect): %d', cat_counts.shape), ...
    sprintf('  Kinetics (rise/decay/plateau): %d', cat_counts.kinetics), ...
    '', ...
    '─── INHERITED ───', ...
    sprintf('  Extreme SNR: %d', cat_counts.extreme_snr), ...
    sprintf('  Correlated: %d', cat_counts.correlation), ...
    sprintf('  Spatial overlap: %d', cat_counts.spatial), ...
    sprintf('  Inactive: %d', cat_counts.inactive), ...
    sprintf('  Baseline drift: %d', cat_counts.baseline)
};

text(0.1, 0.95, summary_text, 'Units', 'normalized', 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'FontName', 'FixedWidth');

sgtitle(sprintf('%s - Refinement Summary (v3.0)', sample_name), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

saveas(fig, fullfile(output_dir, [sample_name '_refinement_summary.png']));
close(fig);

fprintf('    ✓ Generated refinement summary figure\n');

end


function generate_roi_gallery_v3(roi_data, dFF, output_dir, sample_name, frame_rate, config)
% v3.0 ROI Gallery - same as v2.3 but with v3 naming
%
% Layout: [ROI image] [trace] [metrics] per row
% Multiple pages for long recordings

n_rois = size(roi_data.roi_centers, 1);
n_frames = size(dFF, 1);
time_vec = (0:n_frames-1) / frame_rate;
total_duration = time_vec(end);

% Segment parameters
segment_duration = config.gallery_segment_duration;
n_segments = ceil(total_duration / segment_duration);
frames_per_segment = round(segment_duration * frame_rate);

% Page parameters
segments_per_page = 5;
n_pages_per_roi = ceil(n_segments / segments_per_page);

% Get mean image
if isfield(roi_data, 'meanImg')
    meanImg = roi_data.meanImg;
else
    meanImg = zeros(100, 100);
end

% Generate distinct colors
if n_rois <= 20
    roi_colors = lines(n_rois);
else
    roi_colors = hsv(n_rois);
end

% ─────────────────────────────────────────────────────────────────────────
% OVERLAY VIEW
% ─────────────────────────────────────────────────────────────────────────

fig_overlay = figure('Position', [50 50 1400 1000], 'Visible', 'off');

subplot(1, 2, 1);
imagesc(meanImg);
colormap(gca, gray);
axis image; hold on;

if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
    for i = 1:n_rois
        [B, ~] = bwboundaries(roi_data.roi_masks(:,:,i));
        if ~isempty(B)
            plot(B{1}(:,2), B{1}(:,1), 'Color', roi_colors(i,:), 'LineWidth', 2);
        end
        
        center = roi_data.roi_centers(i,:);
        text(center(1), center(2), sprintf('%d', i), 'Color', 'w', 'FontSize', 8, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'BackgroundColor', [roi_colors(i,:) 0.6]);
    end
end

title(sprintf('%s - All %d Refined ROIs', sample_name, n_rois), 'FontSize', 12, 'Interpreter', 'none');
colorbar;

% Legend panel
subplot(1, 2, 2);
axis off; hold on;
xlim([0 1]); ylim([0 1]);

y_pos = 0.95;
text(0.05, y_pos, 'ROI Legend:', 'FontSize', 12, 'FontWeight', 'bold');
y_pos = y_pos - 0.04;

for i = 1:min(n_rois, 20)
    snr_str = '';
    rate_str = '';
    
    if isfield(roi_data, 'snr')
        snr_str = sprintf('SNR=%.1f', roi_data.snr(i));
    end
    if isfield(roi_data, 'event_rate')
        rate_str = sprintf('%.2fHz', roi_data.event_rate(i));
    end
    
    box_x = [0.02, 0.05, 0.05, 0.02];
    box_y = [y_pos-0.012, y_pos-0.012, y_pos+0.012, y_pos+0.012];
    fill(box_x, box_y, roi_colors(i,:), 'EdgeColor', 'k', 'LineWidth', 0.5);
    
    text(0.07, y_pos, sprintf('ROI %d: %s  %s', i, snr_str, rate_str), ...
        'FontSize', 9, 'VerticalAlignment', 'middle');
    
    y_pos = y_pos - 0.045;
end

if n_rois > 20
    text(0.07, y_pos, sprintf('... and %d more ROIs', n_rois - 20), ...
        'FontSize', 9, 'FontAngle', 'italic');
end

sgtitle(sprintf('%s - ROI Overlay Map (v3.0 Refined)', sample_name), 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig_overlay, fullfile(output_dir, sprintf('%s_roi_overlay.png', sample_name)));
close(fig_overlay);

% ─────────────────────────────────────────────────────────────────────────
% GALLERY PAGES
% ─────────────────────────────────────────────────────────────────────────

for roi_num = 1:n_rois
    
    roi_color = roi_colors(roi_num, :);
    full_trace = dFF(:, roi_num);
    
    % Consistent y-limits
    y_min = prctile(full_trace, 0.5);
    y_max = prctile(full_trace, 99.5);
    y_padding = (y_max - y_min) * 0.15;
    if y_padding < 0.01, y_padding = 0.05; end
    y_lim_min = y_min - y_padding;
    y_lim_max = y_max + y_padding;
    
    % ROI image crop
    center = roi_data.roi_centers(roi_num, :);
    crop_size = 50;
    y_range_img = max(1, round(center(2)) - crop_size):min(size(meanImg, 1), round(center(2)) + crop_size);
    x_range_img = max(1, round(center(1)) - crop_size):min(size(meanImg, 2), round(center(1)) + crop_size);
    crop_img = meanImg(y_range_img, x_range_img);
    
    roi_boundary = {};
    if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
        mask_crop = roi_data.roi_masks(y_range_img, x_range_img, roi_num);
        [B, ~] = bwboundaries(mask_crop);
        roi_boundary = B;
    end
    
    % Create pages
    for page = 1:n_pages_per_roi
        
        seg_start_idx = (page - 1) * segments_per_page + 1;
        seg_end_idx = min(page * segments_per_page, n_segments);
        segs_this_page = seg_end_idx - seg_start_idx + 1;
        
        fig_page = figure('Position', [50 50 1800 900], 'Visible', 'off');
        
        for s_idx = 1:segs_this_page
            seg = seg_start_idx + s_idx - 1;
            row = s_idx;
            
            start_frame = (seg - 1) * frames_per_segment + 1;
            end_frame = min(seg * frames_per_segment, n_frames);
            
            if start_frame > n_frames, break; end
            
            seg_frames = start_frame:end_frame;
            seg_trace = dFF(seg_frames, roi_num);
            seg_time = (0:length(seg_frames)-1) / frame_rate;
            
            seg_start_time = (seg - 1) * segment_duration;
            seg_end_time = min(seg * segment_duration, total_duration);
            
            % Panel A: ROI mask crop
            subplot(segments_per_page, 4, (row-1)*4 + 1);
            
            imagesc(crop_img);
            colormap(gca, gray);
            axis image off;
            hold on;
            
            if ~isempty(roi_boundary)
                for k = 1:length(roi_boundary)
                    plot(roi_boundary{k}(:,2), roi_boundary{k}(:,1), ...
                        'Color', roi_color, 'LineWidth', 2.5);
                end
            end
            
            title(sprintf('ROI %d | %d-%ds', roi_num, round(seg_start_time), round(seg_end_time)), ...
                'FontSize', 10, 'Color', roi_color);
            
            % Panel B: Trace
            subplot(segments_per_page, 4, [(row-1)*4 + 2, (row-1)*4 + 3]);
            
            plot(seg_time, seg_trace, 'Color', roi_color, 'LineWidth', 0.8);
            hold on;
            
            baseline = median(full_trace(1:min(100, length(full_trace))));
            yline(baseline, 'k--', 'Alpha', 0.3);
            yline(0, 'k:', 'Alpha', 0.2);
            
            try
                [pks, locs] = findpeaks(seg_trace, 'MinPeakProminence', 0.05, ...
                    'MinPeakDistance', round(0.5 * frame_rate));
                if ~isempty(pks)
                    plot(seg_time(locs), pks, 'rv', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
                end
            catch
            end
            
            xlabel('Time (s)', 'FontSize', 8);
            ylabel('dF/F', 'FontSize', 8);
            xlim([0, segment_duration]);
            ylim([y_lim_min, y_lim_max]);
            set(gca, 'XTick', 0:5:segment_duration);
            grid on;
            set(gca, 'FontSize', 8);
            
            % Panel C: Metrics
            subplot(segments_per_page, 4, (row-1)*4 + 4);
            axis off;
            
            metrics_text = {};
            if isfield(roi_data, 'snr')
                metrics_text{end+1} = sprintf('SNR: %.2f', roi_data.snr(roi_num));
            end
            if isfield(roi_data, 'event_rate')
                metrics_text{end+1} = sprintf('Rate: %.3f Hz', roi_data.event_rate(roi_num));
            end
            if isfield(roi_data, 'event_count')
                metrics_text{end+1} = sprintf('Events: %d', roi_data.event_count(roi_num));
            end
            metrics_text{end+1} = '';
            metrics_text{end+1} = sprintf('Max: %.2f', max(full_trace));
            metrics_text{end+1} = sprintf('Std: %.3f', std(full_trace));
            metrics_text{end+1} = '';
            metrics_text{end+1} = sprintf('Seg Max: %.2f', max(seg_trace));
            
            text(0.1, 0.9, metrics_text, 'Units', 'normalized', 'FontSize', 9, ...
                'VerticalAlignment', 'top');
        end
        
        % Page title
        page_start_time = (seg_start_idx - 1) * segment_duration;
        page_end_time = min(seg_end_idx * segment_duration, total_duration);
        
        if n_pages_per_roi > 1
            sgtitle(sprintf('%s - ROI %d (Page %d/%d: %.0f-%.0fs) [v3.0 Refined]', ...
                sample_name, roi_num, page, n_pages_per_roi, page_start_time, page_end_time), ...
                'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
        else
            sgtitle(sprintf('%s - ROI %d (%.0f-%.0fs) [v3.0 Refined]', ...
                sample_name, roi_num, page_start_time, page_end_time), ...
                'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
        end
        
        % Save
        if n_pages_per_roi > 1
            gallery_path = fullfile(output_dir, sprintf('%s_ROI%02d_gallery_page%d.png', ...
                sample_name, roi_num, page));
        else
            gallery_path = fullfile(output_dir, sprintf('%s_ROI%02d_gallery.png', ...
                sample_name, roi_num));
        end
        
        saveas(fig_page, gallery_path);
        close(fig_page);
    end
end

fprintf('    ✓ Generated overlay + %d ROI galleries\n', n_rois);

end


function generate_refined_analysis_figure_v3(roi_data, output_dir, sample_name, frame_rate)
% v3.0 Analysis Figure - same as v2.3 with v3 naming

n_rois = size(roi_data.roi_centers, 1);
if n_rois == 0
    return;
end

fig = figure('Position', [50 50 1800 1000], 'Visible', 'off');

% Get dFF in correct orientation [time x ROIs]
dFF = double(roi_data.dFF);
if size(dFF, 1) == n_rois && size(dFF, 2) ~= n_rois
    dFF = dFF';
end
n_frames = size(dFF, 1);
time_vec = (0:n_frames-1) / frame_rate;
total_duration = time_vec(end);

% Panel 1: Mean image with refined ROIs
subplot(2, 4, 1);
if isfield(roi_data, 'meanImg')
    imagesc(roi_data.meanImg);
    colormap(gca, gray);
else
    imagesc(zeros(100));
end
axis image; hold on;

if isfield(roi_data, 'snr')
    snr_vals = roi_data.snr(:);
    snr_norm = (snr_vals - min(snr_vals)) / (max(snr_vals) - min(snr_vals) + eps);
    cmap = parula(256);
else
    snr_norm = ones(n_rois, 1) * 0.5;
    cmap = parula(256);
end

if isfield(roi_data, 'roi_masks') && ndims(roi_data.roi_masks) == 3
    for i = 1:n_rois
        color_idx = max(1, min(256, round(snr_norm(i) * 255) + 1));
        roi_color = cmap(color_idx, :);
        
        [B, ~] = bwboundaries(roi_data.roi_masks(:,:,i));
        if ~isempty(B)
            plot(B{1}(:,2), B{1}(:,1), 'Color', roi_color, 'LineWidth', 1.5);
        end
        plot(roi_data.roi_centers(i,1), roi_data.roi_centers(i,2), '+', 'Color', roi_color, 'MarkerSize', 4);
    end
end
title(sprintf('%s - %d ROIs (color=SNR)', sample_name, n_rois), 'Interpreter', 'none', 'FontSize', 10);
cb = colorbar;
cb.Label.String = 'SNR';
if isfield(roi_data, 'snr') && length(snr_vals) > 1 && min(snr_vals) < max(snr_vals)
    caxis([min(snr_vals), max(snr_vals)]);
end

% Panel 2: SNR distribution
subplot(2, 4, 2);
if isfield(roi_data, 'snr')
    histogram(roi_data.snr, 15, 'FaceColor', [0.3 0.7 0.3], 'EdgeColor', 'w');
    xlabel('SNR');
    ylabel('Count');
    xline(mean(roi_data.snr), 'r--', 'LineWidth', 2);
    title(sprintf('SNR (mean=%.1f, med=%.1f)', mean(roi_data.snr), median(roi_data.snr)), 'FontSize', 10);
    grid on;
end

% Panel 3: Event rate distribution
subplot(2, 4, 3);
if isfield(roi_data, 'event_rate')
    histogram(roi_data.event_rate, 15, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
    xlabel('Event Rate (Hz)');
    ylabel('Count');
    xline(mean(roi_data.event_rate), 'r--', 'LineWidth', 2);
    title(sprintf('Event Rate (mean=%.3f Hz)', mean(roi_data.event_rate)), 'FontSize', 10);
    grid on;
end

% Panel 4: Correlation
subplot(2, 4, 4);
if n_rois <= 50
    corr_matrix = corrcoef(dFF);
    imagesc(corr_matrix);
    colormap(gca, create_diverging_colormap());
    caxis([-1, 1]);
    axis square;
    xlabel('ROI');
    ylabel('ROI');
    title(sprintf('Correlation (sync=%.2f)', roi_data.synchrony_score), 'FontSize', 10);
    colorbar;
else
    corr_matrix = corrcoef(dFF);
    upper_tri = triu(true(n_rois), 1);
    corr_vals = corr_matrix(upper_tri);
    
    histogram(corr_vals, 30, 'FaceColor', [0.6 0.4 0.8], 'EdgeColor', 'w');
    xlabel('Pairwise Correlation');
    ylabel('Count');
    xline(mean(corr_vals), 'r--', 'LineWidth', 2);
    title(sprintf('Correlations (sync=%.2f)', roi_data.synchrony_score), 'FontSize', 10);
    xlim([-1, 1]);
    grid on;
end

% Panel 5: Z-score heatmap
subplot(2, 4, 5);

dFF_zscore = zeros(size(dFF));
for i = 1:n_rois
    trace = dFF(:, i);
    trace_mean = mean(trace);
    trace_std = std(trace);
    if trace_std > 0
        dFF_zscore(:, i) = (trace - trace_mean) / trace_std;
    else
        dFF_zscore(:, i) = trace - trace_mean;
    end
end

mean_zscore = mean(dFF_zscore, 1);
[~, sort_idx] = sort(mean_zscore, 'descend');
dFF_zscore_sorted = dFF_zscore(:, sort_idx);

imagesc(time_vec, 1:n_rois, dFF_zscore_sorted');
colormap(gca, create_diverging_colormap());

xlabel('Time (s)');
ylabel('ROI (sorted by activity)');
title('Activity Heatmap (Z-scored)', 'FontSize', 10);
cb = colorbar;
cb.Label.String = 'Z-score';

if total_duration > 200
    set(gca, 'XTick', 0:60:total_duration);
elseif total_duration > 60
    set(gca, 'XTick', 0:30:total_duration);
end

z_max = prctile(abs(dFF_zscore_sorted(:)), 98);
z_max = max(z_max, 1);
caxis([-z_max, z_max]);

% Panel 6: Raw heatmap
subplot(2, 4, 6);

dFF_sorted = dFF(:, sort_idx);
imagesc(time_vec, 1:n_rois, dFF_sorted');
colormap(gca, parula);

xlabel('Time (s)');
ylabel('ROI (sorted by activity)');
title('Activity Heatmap (raw dF/F)', 'FontSize', 10);
cb = colorbar;
cb.Label.String = 'dF/F';

if total_duration > 200
    set(gca, 'XTick', 0:60:total_duration);
elseif total_duration > 60
    set(gca, 'XTick', 0:30:total_duration);
end

dff_p1 = prctile(dFF_sorted(:), 1);
dff_p99 = prctile(dFF_sorted(:), 99);
dff_min = max(dff_p1, -0.2);
dff_max = max(dff_p99, dff_p1 + 0.5);
caxis([dff_min, dff_max]);

% Panel 7: Best traces by SNR
subplot(2, 4, 7);
if isfield(roi_data, 'snr')
    [~, snr_order] = sort(roi_data.snr, 'descend');
else
    snr_order = 1:n_rois;
end
n_show = min(6, n_rois);
show_idx = snr_order(1:n_show);

offset = 0;
colors = lines(n_show);
legend_entries = cell(n_show, 1);

for i = 1:n_show
    trace = dFF(:, show_idx(i));
    plot(time_vec, trace + offset, 'Color', colors(i,:), 'LineWidth', 0.8);
    hold on;
    
    if isfield(roi_data, 'snr')
        legend_entries{i} = sprintf('ROI %d (SNR=%.1f)', show_idx(i), roi_data.snr(show_idx(i)));
    else
        legend_entries{i} = sprintf('ROI %d', show_idx(i));
    end
    
    offset = offset + max(prctile(trace, 99) - prctile(trace, 1), 0.2) * 1.5;
end
xlabel('Time (s)');
ylabel('dF/F (offset)');
title('Best ROIs by SNR (full trace)', 'FontSize', 10);
xlim([0, time_vec(end)]);

if total_duration > 200
    set(gca, 'XTick', 0:60:total_duration);
elseif total_duration > 60
    set(gca, 'XTick', 0:30:total_duration);
end

legend(legend_entries, 'Location', 'eastoutside', 'FontSize', 7);
grid on;

% Panel 8: Summary
subplot(2, 4, 8);
axis off;

sync_status = 'Good';
sync_color = [0 0.6 0];
if roi_data.synchrony_score > 0.5
    sync_status = 'HIGH - check for artifacts';
    sync_color = [0.8 0 0];
elseif roi_data.synchrony_score > 0.3
    sync_status = 'Moderate';
    sync_color = [0.8 0.6 0];
end

stats_text = {
    'REFINED DATA SUMMARY (v3.0)', ...
    '═══════════════════════════', ...
    '', ...
    sprintf('ROIs: %d', n_rois), ...
    sprintf('Duration: %.1f s', n_frames / frame_rate), ...
    '', ...
    sprintf('Mean SNR: %.2f ± %.2f', mean(roi_data.snr), std(roi_data.snr)), ...
    sprintf('Mean Event Rate: %.3f Hz', mean(roi_data.event_rate)), ...
    '', ...
    sprintf('Synchrony: %.2f (%s)', roi_data.synchrony_score, sync_status), ...
    '', ...
    '───────────────────────────', ...
    'v3.0 FEATURES:', ...
    '  Edge contact filter', ...
    '  Shape validation', ...
    '  Kinetics validation', ...
    '  Negative transient check'
};

text(0.05, 0.95, stats_text, 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'FontName', 'FixedWidth');

sgtitle(sprintf('%s - Refined Analysis (v3.0)', sample_name), ...
    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

% Save
saveas(fig, fullfile(output_dir, [sample_name '_refined_analysis.png']));

try
    print(fig, fullfile(output_dir, [sample_name '_refined_analysis.pdf']), '-dpdf', '-bestfit');
catch
end

close(fig);

fprintf('    ✓ Generated refined analysis figure\n');

end


function cmap = create_diverging_colormap(n_colors)
% Create blue-white-red diverging colormap

if nargin < 1
    n_colors = 256;
end

half = floor(n_colors / 2);

blue_to_white = [linspace(0.2, 1, half)', ...
    linspace(0.4, 1, half)', ...
    linspace(0.8, 1, half)'];

white_to_red = [linspace(1, 0.8, half)', ...
    linspace(1, 0.2, half)', ...
    linspace(1, 0.2, half)'];

if mod(n_colors, 2) == 0
    cmap = [blue_to_white; white_to_red];
else
    cmap = [blue_to_white; [1 1 1]; white_to_red];
end

end


function generate_readme_v3(output_dir, summary_table, refine_config, conditions)
% Generate README for v3.0

readme_path = fullfile(output_dir, 'README.txt');
fid = fopen(readme_path, 'w');

fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n');
fprintf(fid, 'ROI REFINEMENT OUTPUT v3.0 - README\n');
fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n\n');

fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));

fprintf(fid, 'v3.0 NEW FEATURES (February 2026)\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, '• EDGE CONTACT FILTER:\n');
fprintf(fid, '  - Removes ROIs with >%.0f%% of perimeter on brain mask edge\n', refine_config.max_edge_contact_fraction * 100);
fprintf(fid, '  - Catches "edge-hugging" artifacts\n\n');
fprintf(fid, '• SHAPE VALIDATION:\n');
fprintf(fid, '  - Compactness > %.2f (removes tails, branches)\n', refine_config.min_compactness);
fprintf(fid, '  - Extent %.2f-%.2f (removes stringy/rectangular)\n', refine_config.min_extent, refine_config.max_extent);
fprintf(fid, '  - Aspect ratio < %.1f (removes very elongated)\n', refine_config.max_aspect_ratio);
fprintf(fid, '  - Solidity > %.2f, Eccentricity < %.2f\n\n', refine_config.min_solidity, refine_config.max_eccentricity);
fprintf(fid, '• KINETICS VALIDATION:\n');
fprintf(fid, '  - Rise/decay ratio < %.1f (GCaMP8f has fast rise, slow decay)\n', refine_config.max_rise_decay_ratio);
fprintf(fid, '  - Decay tau %.2f-%.1f s\n', refine_config.min_decay_tau_s, refine_config.max_decay_tau_s);
fprintf(fid, '  - Max plateau duration: %.1f s\n', refine_config.max_plateau_duration_s);
fprintf(fid, '  - Step artifact detection (instantaneous jumps)\n\n');
fprintf(fid, '• NEGATIVE TRANSIENT CHECK:\n');
fprintf(fid, '  - Removes if z-score < %.1f (dips below baseline)\n\n', refine_config.max_negative_zscore);
fprintf(fid, '• REJECTION EXAMPLES FIGURE:\n');
fprintf(fid, '  - See rejected_examples/ folder\n');
fprintf(fid, '  - Shows WHY each ROI was removed\n\n');

fprintf(fid, 'DATA LOCATIONS\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
for c = 1:length(conditions)
    fprintf(fid, '%s:\n', conditions(c).name);
    fprintf(fid, '  %s\n', conditions(c).output_dir);
    fprintf(fid, '    ├── data/              Refined .mat files\n');
    fprintf(fid, '    ├── figures/           Visualization PNGs\n');
    fprintf(fid, '    ├── logs/              Refinement logs\n');
    fprintf(fid, '    └── rejected_examples/ Why ROIs were removed\n\n');
end

fprintf(fid, 'ALL REFINEMENT CRITERIA\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'v3.0 NEW:\n');
fprintf(fid, '  1. Edge contact > %.0f%% → REMOVE\n', refine_config.max_edge_contact_fraction * 100);
fprintf(fid, '  2. Compactness < %.2f → REMOVE\n', refine_config.min_compactness);
fprintf(fid, '  3. Extent < %.2f or > %.2f → REMOVE\n', refine_config.min_extent, refine_config.max_extent);
fprintf(fid, '  4. Aspect ratio > %.1f → REMOVE\n', refine_config.max_aspect_ratio);
fprintf(fid, '  5. Negative z-score < %.1f → REMOVE\n', refine_config.max_negative_zscore);
fprintf(fid, '  6. Rise/decay ratio > %.1f → REMOVE\n', refine_config.max_rise_decay_ratio);
fprintf(fid, '  7. Step artifacts (≥2 jumps) → REMOVE\n');
fprintf(fid, '  8. Plateau > %.1f s → REMOVE\n\n', refine_config.max_plateau_duration_s);
fprintf(fid, 'INHERITED:\n');
fprintf(fid, '  9. SNR > %.0f → REMOVE\n', refine_config.max_snr);
fprintf(fid, ' 10. Pairwise correlation > %.2f → Remove lower-SNR\n', refine_config.max_pairwise_corr);
fprintf(fid, ' 11. Center distance < %d px → Remove lower-SNR\n', refine_config.min_center_distance_px);
fprintf(fid, ' 12. Solidity < %.2f → REMOVE\n', refine_config.min_solidity);
fprintf(fid, ' 13. Eccentricity > %.2f → REMOVE\n', refine_config.max_eccentricity);
fprintf(fid, ' 14. Event count < %d → REMOVE\n', refine_config.min_event_count);
fprintf(fid, ' 15. Baseline drift > %.2f dF/F → REMOVE\n\n', refine_config.baseline_return_threshold);

fprintf(fid, 'SUMMARY\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
cond_names = unique(summary_table.condition);
for i = 1:length(cond_names)
    cond = cond_names{i};
    rows = strcmp(summary_table.condition, cond);
    
    fprintf(fid, '%s:\n', cond);
    fprintf(fid, '  Files: %d\n', sum(rows));
    fprintf(fid, '  Total ROIs: %d → %d (%.1f%% removed)\n', ...
        sum(summary_table.n_original(rows)), ...
        sum(summary_table.n_refined(rows)), ...
        mean(summary_table.pct_removed(rows)));
    fprintf(fid, '  v3.0 removals:\n');
    fprintf(fid, '    Edge contact: %d\n', sum(summary_table.removed_edge(rows)));
    fprintf(fid, '    Shape: %d\n', sum(summary_table.removed_shape(rows)));
    fprintf(fid, '    Kinetics: %d\n', sum(summary_table.removed_kinetics(rows)));
    fprintf(fid, '    Negative: %d\n', sum(summary_table.removed_negative(rows)));
    fprintf(fid, '    Step: %d\n', sum(summary_table.removed_step(rows)));
    fprintf(fid, '    Plateau: %d\n', sum(summary_table.removed_plateau(rows)));
    fprintf(fid, '  Mean SNR: %.2f\n', mean(summary_table.mean_snr(rows), 'omitnan'));
    fprintf(fid, '  Mean Event Rate: %.3f Hz\n\n', mean(summary_table.mean_event_rate(rows), 'omitnan'));
end

fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n');

fclose(fid);

end