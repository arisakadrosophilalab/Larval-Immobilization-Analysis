function COMPLETE_PIPELINE_RUN_ME_5(varargin)
% COMPLETE_PIPELINE_RUN_ME_v4_5 - Full automated pipeline for larvae calcium imaging
%
% USAGE:
%   COMPLETE_PIPELINE_RUN_ME_v4_5()           % Uses defaults + folder picker
%   COMPLETE_PIPELINE_RUN_ME_v4_5(config)     % Uses custom config struct
%
%   config = pipeline_config();               % Get default config
%   config.conditions(1).name = 'Control';
%   config.conditions(1).tifDir = '/path/to/tif';
%   config.conditions(1).aviDir = '/path/to/avi';
%   COMPLETE_PIPELINE_RUN_ME_v4_5(config);    % Use custom config
%
% INPUTS:
%   config - (optional) Configuration struct from pipeline_config().
%            If omitted, defaults are used. If conditions are empty or
%            paths don't exist, a folder picker opens.
%
% PUBLICATION-READY VERSION - February 2026
%
% ═══════════════════════════════════════════════════════════════════════════
% VERSION 4.5 UPDATES (February 2026)
% ═══════════════════════════════════════════════════════════════════════════
%
% GLOBAL SIGNAL REGRESSION FOR WIDEFIELD IMAGING:
%   Widefield microscopy lacks optical sectioning, causing shared "neuropil"
%   signal to contaminate all ROIs. v4.5 adds optional global signal regression
%   to remove this shared signal and reveal independent neuron activity.
%
%   v4.5 Solutions:
%   1. NEW: Global signal regression (removes shared background fluctuations)
%   2. NEW: Quality metrics stored in roi_data (global_correlation, artifact flags)
%   3. NEW: dFF_raw preserved (before regression) alongside dFF (after)
%   4. NEW: Artifact detection during processing
%
% INHERITED FROM v4.4 (February 2026):
%   - Corrected neuron sizes for Drosophila L2 (gSig=3, max_area=120)
%   - Enhanced figure visualization
%   - Multi-level seed detection
%
% ENHANCED FIGURE VISUALIZATION:
%   1. FILLED ROI OUTLINES: Semi-transparent colored regions
%   2. THICK CONTOURS: 2.5pt lines with white outline for visibility
%   3. ROI NUMBERING: Labels directly on the image
%   4. ZOOMED PANELS: Individual ROI thumbnails with info
%   5. IMPROVED COLORMAP: Better contrast for ROI visualization
%   6. DEDICATED ROI MAP: Clean figure showing just ROI locations
%
% INHERITED FROM v4.3:
%   - Seed-relative thresholds
%   - Non-maximum suppression (NMS)
%   - Savitzky-Golay smoothing
%   - Synchrony detection
%   - Shape metrics (Eccentricity, Solidity)
%   - Y-axis motion correction fix
%
% OUTPUT FILES:
%   - *_rois.mat                  : Full ROI data structure (with quality metrics)
%   - *_roi_metrics_v45.csv       : Per-ROI metrics with shape + activity
%   - *_analysis_v45.png/pdf      : 8-panel diagnostic figure
%   - *_roi_map_v45.png/pdf       : Clean ROI map with numbered outlines
%   - *_roi_gallery_v45.png/pdf   : Zoomed individual ROI thumbnails
%   - *_traces_v45.png/pdf        : Calcium traces (paginated if >16 ROIs)
%   - *_diagnostic_v45.csv        : Diagnostic info when 0 ROIs detected
%   - column_dictionary_v45.txt   : Documentation of all CSV columns
%   - pipeline_summary_v45.csv    : Summary across all files
%
% CONFIGURATION:
%   All parameters (pixel size, frame rate, ROI detection, event detection,
%   QC gating, figure options) are set via the config struct. See
%   pipeline_config() for the full list and documentation.
%
% DEPENDENCIES:
%   - MATLAB R2020b or later (recommended: R2024b)
%   - Bio-Formats (https://www.openmicroscopy.org/bio-formats/)
%   - Image Processing Toolbox
%   - Statistics and Machine Learning Toolbox
%
% For Drosophila larvae L2 neurons with GCaMP8f
% Neuron soma diameter: 2-3 µm (properly sized since v4.4)
% Date: February 2026
%
% See also: pipeline_config

%% =====================================================================
%  CONFIGURATION
%  =====================================================================
if nargin >= 1 && isstruct(varargin{1})
    config = merge_with_defaults(varargin{1});
elseif nargin >= 1
    error(['Unexpected argument (type: %s). Pass a config struct from\n' ...
           'pipeline_config(), or call with no arguments for folder picker.'], ...
           class(varargin{1}));
else
    config = pipeline_config();
end

% If no conditions defined, prompt for folders
if isempty(config.conditions) || ~isfield(config.conditions, 'name') || isempty(config.conditions(1).name)
    config.conditions = prompt_for_conditions();
end

%% =====================================================================
%  PIPELINE START
%  =====================================================================
rng(config.rng_seed, 'twister');
versions = detect_versions();

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║     COMPLETE LARVAE CALCIUM IMAGING PIPELINE v4.5                ║\n');
fprintf('║     GLOBAL SIGNAL REGRESSION FOR WIDEFIELD IMAGING               ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('MATLAB: %s\n', versions.MATLAB_release);
fprintf('Indicator: %s\n', config.indicator);
fprintf('Frame rate: %.4f Hz\n', config.frameRate);
fprintf('Pixel size: %.4f µm\n', config.pixelSize_um);
fprintf('\n');
fprintf('═══ NEURON SIZE PARAMETERS (from v4.4) ═══\n');
fprintf('  gSig = %d px (%.2f µm radius, %.2f µm diameter)\n', ...
    config.roi.gSig, config.roi.gSig * config.pixelSize_um, 2 * config.roi.gSig * config.pixelSize_um);
fprintf('  min_area = %d px (%.2f µm², ~%.2f µm diameter)\n', ...
    config.roi.min_area_px, config.roi.min_area_px * config.pixelSize_um^2, ...
    2 * sqrt(config.roi.min_area_px * config.pixelSize_um^2 / pi));
fprintf('  max_area = %d px (%.2f µm², ~%.2f µm diameter)\n', ...
    config.roi.max_area_px, config.roi.max_area_px * config.pixelSize_um^2, ...
    2 * sqrt(config.roi.max_area_px * config.pixelSize_um^2 / pi));
target_area = pi * (config.roi.gSig * config.roi.target_area_factor)^2;
fprintf('  target_area = %.0f px (%.2f µm², ~%.2f µm diameter)\n', ...
    target_area, target_area * config.pixelSize_um^2, ...
    2 * sqrt(target_area * config.pixelSize_um^2 / pi));
fprintf('═══════════════════════════════════════════════\n');
fprintf('\n');
fprintf('═══ v4.5 GLOBAL SIGNAL REGRESSION ═══\n');
fprintf('  apply_global_regression = %s\n', mat2str(config.roi.apply_global_regression));
fprintf('  method = %s\n', config.roi.global_regression_method);
fprintf('  Removes shared neuropil signal in widefield imaging\n');
fprintf('═══════════════════════════════════════════════\n');
fprintf('\n');

% v4.4: Validate parameters before processing
validate_config(config);

setup_bioformats();

all_results = table();
all_thresholds = {};

for condIdx = 1:length(config.conditions)
    cond = config.conditions(condIdx);
    
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  CONDITION: %-52s ║\n', cond.name);
    fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
    fprintf('\n');
    fprintf('TIF folder: %s\n', cond.tifDir);
    fprintf('AVI folder: %s\n', cond.aviDir);
    fprintf('\n');
    
    if ~isfolder(cond.tifDir)
        fprintf('[SKIP] TIF directory not found: %s\n', cond.tifDir);
        continue;
    end
    
    [eligible_files, qc_info] = load_qc_eligible_files(cond, config);
    
    if isempty(eligible_files)
        if config.require_qc_results
            fprintf('[SKIP] No eligible files found\n');
            continue;
        else
            tifFiles = find_input_files(cond.tifDir, config.file_pattern);
            if isempty(tifFiles)
                fprintf('[SKIP] No S*.tif files found\n');
                continue;
            end
            eligible_files = {tifFiles.name};
        end
    end
    
    fprintf('Files to process: %d\n', length(eligible_files));
    if ~isempty(qc_info) && ~isempty(qc_info.source)
        fprintf('QC source: %s\n', qc_info.source);
    end
    fprintf('\n');
    
    outputDirs = create_output_dirs(cond.tifDir);
    
    % Generate column dictionary at start of each condition
    write_column_dictionary(outputDirs.rois, config);
    
    logFile = fullfile(cond.tifDir, sprintf('pipeline_log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS')));
    logFid = fopen(logFile, 'w');
    if logFid < 0, logFid = -1; end
    
    logf(logFid, '════════════════════════════════════════════════════════════════\n');
    logf(logFid, 'ROI Detection Pipeline v4.5 Started: %s\n', datestr(now));
    logf(logFid, 'Condition: %s\n', cond.name);
    logf(logFid, 'Files to process: %d\n', length(eligible_files));
    logf(logFid, 'Detection params: gSig=%d, min_area=%d, max_area=%d\n', ...
        config.roi.gSig, config.roi.min_area_px, config.roi.max_area_px);
    logf(logFid, '════════════════════════════════════════════════════════════════\n\n');
    
    condition_thresholds = {};
    
    for idx = 1:length(eligible_files)
        basename = regexprep(eligible_files{idx}, '\.(tif|tiff|avi)$', '', 'ignorecase');
        
        tifFile = find_tif_file(cond.tifDir, basename);
        if isempty(tifFile)
            fprintf('[SKIP] TIF file not found for: %s\n', basename);
            continue;
        end
        
        fprintf('══════════ FILE %d/%d: %s ══════════\n', idx, length(eligible_files), basename);
        logf(logFid, '\n─── %s ───\n', basename);
        
        shiftsCSV = fullfile(cond.aviDir, [basename '_shifts.csv']);
        if ~exist(shiftsCSV, 'file')
            fprintf('  [SKIP] No shifts file found\n\n');
            continue;
        end
        
        try
            [result, thresh_info] = process_single_file(tifFile, basename, config, outputDirs, cond, versions, logFid);
            all_results = append_result(all_results, result);
            condition_thresholds{end+1} = thresh_info;
            fprintf('✓ Complete: %d ROIs (SNR: %.2f, Rate: %.3f Hz)\n\n', ...
                result.nROIs, result.meanSNR, result.meanEventRate);
        catch ME
            fprintf('✗ Error: %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('  at %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
            fprintf('\n');
            logf(logFid, '  ✗ Error: %s\n', ME.message);
        end
    end
    
    % Write thresholds TSV for this condition
    if ~isempty(condition_thresholds)
        write_thresholds_tsv(outputDirs.rois, condition_thresholds, cond.name);
        all_thresholds = [all_thresholds, condition_thresholds];
    end
    
    if ~isempty(all_results)
        cond_results = all_results(strcmp(all_results.condition, cond.name), :);
        if height(cond_results) > 0
            write_condition_summary(cond_results, cond.tifDir, cond.name, config, versions);
        end
    end
    
    if logFid > 0, fclose(logFid); end
end

if ~isempty(all_results)
    generate_final_summary(all_results, config.conditions(1).tifDir, config, versions);
end

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                      PIPELINE COMPLETE                           ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
if ~isempty(all_results)
    fprintf('Total files processed: %d\n', height(all_results));
    fprintf('Total ROIs detected: %d\n', sum(all_results.nROIs));
end
fprintf('\n');

end

%% =====================================================================
%  WRITE COLUMN DICTIONARY
%  =====================================================================
function write_column_dictionary(outputDir, config)
% Generate documentation file explaining all CSV columns

dictFile = fullfile(outputDir, 'column_dictionary_v45.txt');
fid = fopen(dictFile, 'w');
if fid < 0
    warning('Could not create column dictionary file');
    return;
end

fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n');
fprintf(fid, 'COLUMN DICTIONARY - ROI Detection Pipeline v4.5\n');
fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '═══════════════════════════════════════════════════════════════════════════════\n\n');

fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'NEURON SIZE CORRECTIONS (v4.4+)\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'Pixel size: %.4f µm\n', config.pixelSize_um);
fprintf(fid, 'gSig = %d pixels (%.2f µm radius)\n', config.roi.gSig, config.roi.gSig * config.pixelSize_um);
fprintf(fid, 'min_area = %d pixels (%.2f µm² → ~%.2f µm diameter)\n', ...
    config.roi.min_area_px, config.roi.min_area_px * config.pixelSize_um^2, ...
    2 * sqrt(config.roi.min_area_px * config.pixelSize_um^2 / pi));
fprintf(fid, 'max_area = %d pixels (%.2f µm² → ~%.2f µm diameter)\n', ...
    config.roi.max_area_px, config.roi.max_area_px * config.pixelSize_um^2, ...
    2 * sqrt(config.roi.max_area_px * config.pixelSize_um^2 / pi));
fprintf(fid, 'Target: Drosophila L2 neuron somata (2-3 µm diameter)\n\n');

fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'FILE: *_roi_metrics_v45.csv\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n\n');

fprintf(fid, 'COLUMN              UNIT        DESCRIPTION\n');
fprintf(fid, '─────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'ROI_ID              Integer     Unique identifier for each ROI (1-indexed)\n');
fprintf(fid, 'Center_X            Pixels      X coordinate of ROI centroid (column)\n');
fprintf(fid, 'Center_Y            Pixels      Y coordinate of ROI centroid (row)\n');
fprintf(fid, 'Area_px             Pixels²     Total area of ROI mask\n');
fprintf(fid, 'Area_um2            µm²         Total area in physical units\n');
fprintf(fid, 'Diameter_um         µm          Equivalent circular diameter\n');
fprintf(fid, 'Eccentricity        0-1         Shape: 0=circular, 1=line-like\n');
fprintf(fid, 'Solidity            0-1         Shape: 1=convex, <1=irregular\n');
fprintf(fid, 'Detection_Level     Integer     Detection method (1-6)\n');
fprintf(fid, 'Detection_Name      String      Human-readable detection level name\n');
fprintf(fid, 'SNR                 Ratio       Signal-to-noise ratio\n');
fprintf(fid, 'Event_Rate_Hz       Hz          Detected calcium events per second\n');
fprintf(fid, 'Event_Amplitude     dF/F        Mean amplitude of detected events\n');
fprintf(fid, 'Event_Count         Integer     Total number of detected events\n');
fprintf(fid, 'Mean_dFF            dF/F        Mean of dF/F trace\n');
fprintf(fid, 'Max_dFF             dF/F        Maximum dF/F value\n');
fprintf(fid, 'Std_dFF             dF/F        Standard deviation of dF/F\n');
fprintf(fid, 'P25_dFF             dF/F        25th percentile\n');
fprintf(fid, 'P50_dFF             dF/F        50th percentile (median)\n');
fprintf(fid, 'P75_dFF             dF/F        75th percentile\n');
fprintf(fid, 'P90_dFF             dF/F        90th percentile\n');
fprintf(fid, 'P95_dFF             dF/F        95th percentile\n\n');

fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'MOTION CORRECTION CONVENTION - DO NOT CHANGE\n');
fprintf(fid, '───────────────────────────────────────────────────────────────────────────────\n');
fprintf(fid, 'Correction applied as: imtranslate(img, [-dx, +dy])\n\n');
fprintf(fid, 'Why the signs are ASYMMETRIC (and correct):\n');
fprintf(fid, '  • Tracker uses Cartesian coordinates: dx>0 = moved RIGHT, dy>0 = moved UP\n');
fprintf(fid, '  • imtranslate uses image coordinates: +tx = shift RIGHT, +ty = shift DOWN\n');
fprintf(fid, '  • To correct motion, shift OPPOSITE to movement direction:\n');
fprintf(fid, '      X: sample moved right → shift left → use -dx ✓\n');
fprintf(fid, '      Y: sample moved up → shift down → use +dy (NOT -dy!) ✓\n');
fprintf(fid, '  • The +dy works because imtranslate +ty = DOWN, which is what we want\n');
fprintf(fid, '    when dy>0 (sample moved UP)\n\n');
fprintf(fid, 'WARNING: Do NOT "fix" this to [-dx, -dy] - that would BREAK Y-axis correction!\n\n');

fprintf(fid, 'Pipeline version: 4.5\n');
fprintf(fid, 'Date: February 2026\n');

fclose(fid);
fprintf('    ✓ Column dictionary: %s\n', dictFile);

end

%% =====================================================================
%  WRITE THRESHOLDS TSV
%  =====================================================================
function write_thresholds_tsv(outputDir, thresh_data, condName)

tsvFile = fullfile(outputDir, sprintf('roi_thresholds_%s_v45.tsv', lower(strrep(condName, ' ', '_'))));
fid = fopen(tsvFile, 'w');
if fid < 0
    warning('Could not create thresholds TSV file');
    return;
end

fprintf(fid, 'Filename\tCondition\tnROIs\tCn_median\tCn_range\tCn_max\tPNR_median\tPNR_max\t');
fprintf(fid, 'Cn_thresh_used\tPNR_thresh_used\tDetection_method\tCn_uniform\tMotion_quality\n');

for i = 1:length(thresh_data)
    t = thresh_data{i};
    fprintf(fid, '%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.3f\t%.2f\t%s\t%s\t%s\n', ...
        t.filename, t.condition, t.nROIs, ...
        t.cn_median, t.cn_range, t.cn_max, ...
        t.pnr_median, t.pnr_max, ...
        t.cn_thresh_used, t.pnr_thresh_used, ...
        t.detection_method, t.cn_uniform, t.motion_quality);
end

fclose(fid);
fprintf('    ✓ Thresholds TSV: %s\n', tsvFile);

end

%% =====================================================================
%  LOAD QC-ELIGIBLE FILES
%  =====================================================================
function [eligible_files, qc_info] = load_qc_eligible_files(cond, config)

eligible_files = {};
qc_info = struct('source', '', 'column_used', '', 'total_files', 0, 'eligible_files', 0);

if ~isfolder(cond.aviDir)
    return;
end

% Option 1: roi_eligible_list.txt
if strcmpi(config.qc_gate_column, 'ROI_ELIGIBLE')
    listFile = fullfile(cond.aviDir, 'roi_eligible_list.txt');
    if exist(listFile, 'file')
        fid = fopen(listFile, 'r');
        if fid > 0
            stems = textscan(fid, '%s');
            fclose(fid);
            eligible_files = stems{1};
            qc_info.source = listFile;
            qc_info.column_used = 'ROI_ELIGIBLE';
            qc_info.eligible_files = length(eligible_files);
            fprintf('Loaded %d eligible files from roi_eligible_list.txt\n', length(eligible_files));
            return;
        end
    end
end

% Option 2: roi_candidates.csv
candidatesFile = fullfile(cond.aviDir, 'roi_candidates.csv');
if exist(candidatesFile, 'file')
    try
        T = readtable(candidatesFile);
        qc_info.source = candidatesFile;
        qc_info.total_files = height(T);
        
        if ismember(config.qc_gate_column, T.Properties.VariableNames)
            pass_idx = logical(T.(config.qc_gate_column));
            qc_info.column_used = config.qc_gate_column;
        elseif ismember('ROI_ELIGIBLE', T.Properties.VariableNames)
            pass_idx = logical(T.ROI_ELIGIBLE);
            qc_info.column_used = 'ROI_ELIGIBLE';
        else
            pass_idx = true(height(T), 1);
            qc_info.column_used = 'none';
        end
        
        if ismember('Stem', T.Properties.VariableNames)
            eligible_files = T.Stem(pass_idx);
        else
            eligible_files = T{pass_idx, 1};
        end
        
        if ~iscell(eligible_files)
            eligible_files = cellstr(string(eligible_files));
        end
        
        qc_info.eligible_files = length(eligible_files);
        fprintf('Loaded %d/%d eligible files from roi_candidates.csv\n', length(eligible_files), qc_info.total_files);
        return;
    catch
    end
end

end

%% =====================================================================
%  FIND TIF FILE
%  =====================================================================
function tifFile = find_tif_file(tifDir, basename)
tifFile = '';
extensions = {'.tif', '.tiff', '.TIF', '.TIFF'};
for i = 1:length(extensions)
    candidate = fullfile(tifDir, [basename extensions{i}]);
    if exist(candidate, 'file')
        tifFile = candidate;
        return;
    end
end
end

%% =====================================================================
%  MAIN PROCESSING FUNCTION
%  =====================================================================
function [result, thresh_info] = process_single_file(filepath, basename, config, outputDirs, cond, versions, logFid)

result = struct();
result.filename = basename;
result.condition = string(cond.name);
result.nROIs = 0;
result.meanSNR = NaN;
result.meanEventRate = NaN;
result.meanEventAmplitude = NaN;
result.synchronyScore = NaN;  % v4.4: Initialize synchrony score
result.motionCorrected = false;
result.motionQuality = "unknown";
result.qc_gate = config.qc_gate_column;
result.processed_utc = string(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));

thresh_info = struct();
thresh_info.filename = basename;
thresh_info.condition = cond.name;

motionFile = fullfile(outputDirs.motion, [basename '_rigid.tif']);
shiftsFile = fullfile(outputDirs.motion, [basename '_rigid_shifts.mat']);
roiFile = fullfile(outputDirs.rois, [basename '_rois.mat']);
figFile = fullfile(outputDirs.figures, [basename '_analysis_v45.png']);
csvFile = fullfile(outputDirs.rois, [basename '_roi_metrics_v45.csv']);

% Skip if already processed
if exist(roiFile, 'file') && config.skip_existing
    fprintf('  [SKIP] Already processed - loading existing results\n');
    logf(logFid, '  [SKIP] Using existing ROI file\n');
    
    loaded = load(roiFile);
    roi_data = loaded.roi_data;
    
    result.nROIs = size(roi_data.roi_centers, 1);
    result.motionCorrected = true;
    
    if isfield(roi_data, 'provenance') && isfield(roi_data.provenance, 'motion_quality')
        result.motionQuality = string(roi_data.provenance.motion_quality);
    else
        result.motionQuality = "existing";
    end
    
    if result.nROIs > 0
        result.meanSNR = mean(roi_data.snr);
        result.meanEventRate = mean(roi_data.event_rate);
        if isfield(roi_data, 'event_amplitude')
            result.meanEventAmplitude = mean(roi_data.event_amplitude);
        end
        if isfield(roi_data, 'synchrony_score')
            result.synchronyScore = roi_data.synchrony_score;
        else
            result.synchronyScore = NaN;
        end
    else
        result.meanSNR = 0;
        result.meanEventRate = 0;
        result.meanEventAmplitude = 0;
        result.synchronyScore = NaN;
    end
    
    % Populate threshold info
    if isfield(roi_data, 'detection_info')
        info = roi_data.detection_info;
        thresh_info.nROIs = result.nROIs;
        thresh_info.cn_median = info.cn_median;
        thresh_info.cn_range = info.cn_range;
        thresh_info.cn_max = info.cn_max;
        thresh_info.pnr_median = info.pnr_median;
        thresh_info.pnr_max = info.pnr_max;
        
        if isfield(info, 'threshold_levels') && ~isempty(info.threshold_levels)
            thresh_info.cn_thresh_used = info.threshold_levels{end}.cn;
            thresh_info.pnr_thresh_used = info.threshold_levels{end}.pnr;
        else
            thresh_info.cn_thresh_used = config.roi.min_corr;
            thresh_info.pnr_thresh_used = config.roi.min_pnr;
        end
        
        thresh_info.detection_method = 'existing';
        thresh_info.cn_uniform = 'false';
        if isfield(info, 'cn_uniform') && info.cn_uniform
            thresh_info.cn_uniform = 'true';
        end
    else
        thresh_info.nROIs = result.nROIs;
        thresh_info.cn_median = NaN;
        thresh_info.cn_range = NaN;
        thresh_info.cn_max = NaN;
        thresh_info.pnr_median = NaN;
        thresh_info.pnr_max = NaN;
        thresh_info.cn_thresh_used = config.roi.min_corr;
        thresh_info.pnr_thresh_used = config.roi.min_pnr;
        thresh_info.detection_method = 'existing';
        thresh_info.cn_uniform = 'unknown';
    end
    thresh_info.motion_quality = char(result.motionQuality);
    
    % Generate missing figures (check for v44, v43, AND generic versions)
    % v4.4 figures: *_analysis_v45.png, *_roi_map_v45.png, etc.
    % v4.3 figures: *_analysis_v43.png, *_roi_map_v43.png, etc.
    % Generic: *_analysis.png, *_roi_map.png, etc.
    % Set config.force_regenerate_figs = true to regenerate ALL figures
    
    figFile_v44 = fullfile(outputDirs.figures, [basename '_analysis_v45.png']);
    figFile_v43 = fullfile(outputDirs.figures, [basename '_analysis_v43.png']);
    figFile_old = fullfile(outputDirs.figures, [basename '_analysis.png']);
    has_analysis_fig = exist(figFile_v44, 'file') || exist(figFile_v43, 'file') || exist(figFile_old, 'file');
    
    if config.generate_figures && (~has_analysis_fig || config.force_regenerate_figs)
        if config.force_regenerate_figs
            fprintf('  Regenerating analysis figure (force mode)...\n');
        else
            fprintf('  Generating missing analysis figure...\n');
        end
        generate_analysis_figures_v45(roi_data, outputDirs.figures, basename, config);
    end
    
    if config.generate_roi_map && result.nROIs > 0
        roiMapFile_v44 = fullfile(outputDirs.figures, [basename '_roi_map_v45.png']);
        roiMapFile_v43 = fullfile(outputDirs.figures, [basename '_roi_map_v43.png']);
        roiMapFile_old = fullfile(outputDirs.figures, [basename '_roi_map.png']);
        has_roi_map = exist(roiMapFile_v44, 'file') || exist(roiMapFile_v43, 'file') || exist(roiMapFile_old, 'file');
        if ~has_roi_map || config.force_regenerate_figs
            if config.force_regenerate_figs
                fprintf('  Regenerating ROI map (force mode)...\n');
            else
                fprintf('  Generating missing ROI map...\n');
            end
            generate_roi_map_v45(roi_data, outputDirs.figures, basename, config);
        end
    end
    
    if config.generate_roi_gallery && result.nROIs > 0
        % Check for gallery - handle both paginated (_1.png) and non-paginated (.png)
        galleryFile_v44 = fullfile(outputDirs.figures, [basename '_roi_gallery_v45.png']);
        galleryFile_v44_p1 = fullfile(outputDirs.figures, [basename '_roi_gallery_v45_1.png']);
        galleryFile_v43 = fullfile(outputDirs.figures, [basename '_roi_gallery_v43.png']);
        galleryFile_v43_p1 = fullfile(outputDirs.figures, [basename '_roi_gallery_v43_1.png']);
        galleryFile_old = fullfile(outputDirs.figures, [basename '_roi_gallery.png']);
        galleryFile_old_p1 = fullfile(outputDirs.figures, [basename '_roi_gallery_1.png']);
        has_gallery = exist(galleryFile_v44, 'file') || exist(galleryFile_v44_p1, 'file') || ...
                      exist(galleryFile_v43, 'file') || exist(galleryFile_v43_p1, 'file') || ...
                      exist(galleryFile_old, 'file') || exist(galleryFile_old_p1, 'file');
        if ~has_gallery || config.force_regenerate_figs
            if config.force_regenerate_figs
                fprintf('  Regenerating ROI gallery (force mode)...\n');
            else
                fprintf('  Generating missing ROI gallery...\n');
            end
            generate_roi_gallery_v45(roi_data, outputDirs.figures, basename, config);
        end
    end
    
    if config.generate_trace_fig && result.nROIs > 0
        % Check for traces - handle both paginated (_1.png) and non-paginated (.png)
        traceFig_v44 = fullfile(outputDirs.figures, [basename '_traces_v45.png']);
        traceFig_v44_p1 = fullfile(outputDirs.figures, [basename '_traces_v45_1.png']);
        traceFig_v43 = fullfile(outputDirs.figures, [basename '_traces_v43.png']);
        traceFig_v43_p1 = fullfile(outputDirs.figures, [basename '_traces_v43_1.png']);
        traceFig_old = fullfile(outputDirs.figures, [basename '_traces.png']);
        traceFig_old_p1 = fullfile(outputDirs.figures, [basename '_traces_1.png']);
        has_traces = exist(traceFig_v44, 'file') || exist(traceFig_v44_p1, 'file') || ...
                     exist(traceFig_v43, 'file') || exist(traceFig_v43_p1, 'file') || ...
                     exist(traceFig_old, 'file') || exist(traceFig_old_p1, 'file');
        if ~has_traces || config.force_regenerate_figs
            if config.force_regenerate_figs
                fprintf('  Regenerating trace figure (force mode)...\n');
            else
                fprintf('  Generating missing trace figure...\n');
            end
            generate_trace_figures_v45(roi_data, outputDirs.figures, basename, config);
        end
    end
    
    csvFile_v44 = fullfile(outputDirs.rois, [basename '_roi_metrics_v45.csv']);
    csvFile_v43 = fullfile(outputDirs.rois, [basename '_roi_metrics_v43.csv']);
    csvFile_old = fullfile(outputDirs.rois, [basename '_roi_metrics.csv']);
    has_csv = exist(csvFile_v44, 'file') || exist(csvFile_v43, 'file') || exist(csvFile_old, 'file');
    if ~has_csv || config.force_regenerate_figs
        if config.force_regenerate_figs
            fprintf('  Regenerating CSV (force mode)...\n');
        else
            fprintf('  Generating missing CSV...\n');
        end
        export_roi_csv_v45(outputDirs.rois, basename, roi_data, config);
    end
    
    return;
end

% Step 1: Motion Correction
fprintf('  Step 1: Rigid Motion Correction\n');
logf(logFid, '  Step 1: Motion correction\n');

if exist(motionFile, 'file') && config.skip_existing
    fprintf('    ✓ Using existing motion-corrected file\n');
    result.motionCorrected = true;
    
    if exist(shiftsFile, 'file')
        try
            s = load(shiftsFile, 'shifts');
            result.motionQuality = string(s.shifts.motion_quality);
        catch
            result.motionQuality = "existing";
        end
    else
        result.motionQuality = "existing";
    end
else
    shiftsCSV = fullfile(cond.aviDir, [basename '_shifts.csv']);
    fprintf('    Applying rigid shifts...\n');
    
    [success, motionQuality] = apply_rigid_correction_bigtiff(filepath, motionFile, shiftsCSV, shiftsFile, config, logFid);
    result.motionCorrected = success;
    result.motionQuality = motionQuality;
    
    if ~success
        error('Motion correction failed');
    end
end

thresh_info.motion_quality = char(result.motionQuality);

% Step 2: ROI Detection
fprintf('  Step 2: ROI Detection (v4.4 corrected sizes)\n');
logf(logFid, '  Step 2: ROI detection (v4.5)\n');

fprintf('    Detecting ROIs...\n');
roi_data = detect_rois_with_masking(motionFile, config, logFid);

% Populate threshold info
info = roi_data.detection_info;
thresh_info.nROIs = size(roi_data.roi_centers, 1);
thresh_info.cn_median = info.cn_median;
thresh_info.cn_range = info.cn_range;
thresh_info.cn_max = info.cn_max;
thresh_info.pnr_median = info.pnr_median;
thresh_info.pnr_max = info.pnr_max;

if isfield(info, 'threshold_levels') && ~isempty(info.threshold_levels)
    thresh_info.cn_thresh_used = info.threshold_levels{end}.cn;
    thresh_info.pnr_thresh_used = info.threshold_levels{end}.pnr;
else
    thresh_info.cn_thresh_used = config.roi.min_corr;
    thresh_info.pnr_thresh_used = config.roi.min_pnr;
end

if isfield(info, 'cn_uniform') && info.cn_uniform
    thresh_info.detection_method = 'pnr_only';
    thresh_info.cn_uniform = 'true';
elseif thresh_info.nROIs > 0 && isfield(info, 'roi_levels')
    levels = info.roi_levels;
    if any(levels <= 3)
        thresh_info.detection_method = 'cn_pnr';
    elseif any(levels == 4)
        thresh_info.detection_method = 'pnr_only';
    else
        thresh_info.detection_method = 'fallback';
    end
    thresh_info.cn_uniform = 'false';
else
    thresh_info.detection_method = 'none';
    thresh_info.cn_uniform = 'false';
end

roi_data.provenance = struct();
roi_data.provenance.created_utc = result.processed_utc;
roi_data.provenance.matlab_release = versions.MATLAB_release;
roi_data.provenance.pipeline_version = versions.pipeline_version;
roi_data.provenance.config = config.roi;
roi_data.provenance.event_config = config.events;
roi_data.provenance.motion_quality = result.motionQuality;
roi_data.provenance.growing_method = 'seed_relative_v44';

save(roiFile, 'roi_data', '-v7.3');
logf(logFid, '    Saved: %s\n', roiFile);

result.nROIs = size(roi_data.roi_centers, 1);
if result.nROIs > 0
    result.meanSNR = mean(roi_data.snr);
    result.meanEventRate = mean(roi_data.event_rate);
    if isfield(roi_data, 'event_amplitude')
        result.meanEventAmplitude = mean(roi_data.event_amplitude);
    end
    if isfield(roi_data, 'synchrony_score')
        result.synchronyScore = roi_data.synchrony_score;
    else
        result.synchronyScore = NaN;
    end
else
    result.meanSNR = 0;
    result.meanEventRate = 0;
    result.meanEventAmplitude = 0;
    result.synchronyScore = NaN;
end

% Step 3: Figures
if config.generate_figures && (result.nROIs > 0 || config.save_diagnostic_figures)
    fprintf('  Step 3: Generating figures (v4.4 enhanced)\n');
    generate_analysis_figures_v45(roi_data, outputDirs.figures, basename, config);
    
    if config.generate_roi_map && result.nROIs > 0
        generate_roi_map_v45(roi_data, outputDirs.figures, basename, config);
    end
    
    if config.generate_roi_gallery && result.nROIs > 0
        generate_roi_gallery_v45(roi_data, outputDirs.figures, basename, config);
    end
    
    if config.generate_trace_fig && result.nROIs > 0
        generate_trace_figures_v45(roi_data, outputDirs.figures, basename, config);
    end
end

% Step 4: CSV
export_roi_csv_v45(outputDirs.rois, basename, roi_data, config);

end

%% =====================================================================
%  RIGID MOTION CORRECTION
%  =====================================================================
function [success, motionQuality] = apply_rigid_correction_bigtiff(inputFile, outputFile, shiftsCSV, shiftsMatFile, config, logFid)

success = false;
motionQuality = "unknown";

try
    shiftsTable = readtable(shiftsCSV);
    
    if ismember('dx_px', shiftsTable.Properties.VariableNames)
        dx = shiftsTable.dx_px;
        dy = shiftsTable.dy_px;
    else
        error('Shifts CSV missing dx_px/dy_px columns');
    end
    
    fprintf('    Loaded %d shifts: dx=[%.1f, %.1f], dy=[%.1f, %.1f] px\n', ...
        length(dx), min(dx), max(dx), min(dy), max(dy));
    logf(logFid, '    Shifts: dx=[%.1f,%.1f], dy=[%.1f,%.1f] px\n', min(dx), max(dx), min(dy), max(dy));
    
    max_shift = max(max(abs(dx)), max(abs(dy)));
    shift_range = max(max(dx) - min(dx), max(dy) - min(dy));
    
    if max_shift > config.motion.max_reasonable_shift
        fprintf('    ⚠ WARNING: LARGE SHIFTS (max=%.1f px > threshold=%d px)\n', ...
            max_shift, config.motion.max_reasonable_shift);
        motionQuality = "excessive_motion";
    elseif shift_range > config.motion.max_reasonable_shift
        motionQuality = "high_motion";
    else
        motionQuality = "good";
    end
    
    reader = bfGetReader(inputFile);
    reader.setSeries(0);
    width = reader.getSizeX();
    height = reader.getSizeY();
    nFrames = reader.getImageCount();
    
    fprintf('    Input: %dx%d, %d frames\n', height, width, nFrames);
    
    if length(dx) ~= nFrames
        warning('Frame mismatch: shifts=%d, TIFF=%d', length(dx), nFrames);
        if length(dx) < nFrames
            dx = [dx; repmat(dx(end), nFrames - length(dx), 1)];
            dy = [dy; repmat(dy(end), nFrames - length(dy), 1)];
        else
            dx = dx(1:nFrames);
            dy = dy(1:nFrames);
        end
    end
    
    if exist(outputFile, 'file')
        delete(outputFile);
        pause(0.1);
    end
    
    chunkSize = config.motion.chunk_size;
    nChunks = ceil(nFrames / chunkSize);
    
    for chunk = 1:nChunks
        startIdx = (chunk - 1) * chunkSize + 1;
        endIdx = min(chunk * chunkSize, nFrames);
        chunkFrames = endIdx - startIdx + 1;
        
        fprintf('    Chunk %d/%d (frames %d-%d)...\n', chunk, nChunks, startIdx, endIdx);
        
        corrected = zeros(height, width, chunkFrames, 'uint16');
        
        for i = 1:chunkFrames
            frameIdx = startIdx + i - 1;
            img = single(bfGetPlane(reader, frameIdx));
            
            % Motion correction: imtranslate(img, [-dx, +dy])
            img_corrected = imtranslate(img, [-dx(frameIdx), +dy(frameIdx)], 'FillValues', NaN, 'OutputView', 'same');
            
            nan_mask = isnan(img_corrected);
            if any(nan_mask(:))
                img_corrected(nan_mask) = nanmean(img_corrected(:));
            end
            corrected(:,:,i) = uint16(img_corrected);
        end
        
        write_bigtiff_chunk(outputFile, corrected, chunk == 1);
        clear corrected;
    end
    
    reader.close();
    
    if config.save_shifts
        shifts = struct('dx', dx, 'dy', dy, 'source', shiftsCSV, ...
            'max_shift', max_shift, 'motion_quality', motionQuality, ...
            'correction_convention', 'imtranslate(img, [-dx, +dy])');
        save(shiftsMatFile, 'shifts', '-v7.3');
    end
    
    fprintf('    ✓ Motion correction complete (quality: %s)\n', motionQuality);
    success = true;
    
catch ME
    fprintf('    ✗ Error: %s\n', ME.message);
    logf(logFid, '    ✗ Error: %s\n', ME.message);
    if exist('reader', 'var')
        try reader.close(); catch; end
    end
end

end

%% =====================================================================
%  BIGTIFF WRITING
%  =====================================================================
function write_bigtiff_chunk(filename, data, is_first_chunk)

max_retries = 5;
base_pause = 0.2;

for i = 1:size(data, 3)
    frame = uint16(data(:,:,i));
    
    if is_first_chunk && i == 1
        t = Tiff(filename, 'w8');
        t.setTag('ImageWidth', size(frame, 2));
        t.setTag('ImageLength', size(frame, 1));
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('RowsPerStrip', size(frame, 1));
        t.setTag('Software', 'MATLAB Calcium Imaging Pipeline v4.5');
        t.write(frame);
        t.close();
        pause(0.05);
    else
        success = false;
        for attempt = 1:max_retries
            try
                t = Tiff(filename, 'r+');
                success = true;
                break;
            catch
                pause(base_pause * (2 ^ (attempt - 1)));
            end
        end
        
        if ~success
            pause(0.5);
            t = Tiff(filename, 'r+');
        end
        
        while ~t.lastDirectory(), t.nextDirectory(); end
        t.writeDirectory();
        t.setTag('ImageWidth', size(frame, 2));
        t.setTag('ImageLength', size(frame, 1));
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('RowsPerStrip', size(frame, 1));
        t.write(frame);
        t.close();
        
        if mod(i, 25) == 0, pause(0.02); end
    end
end
end

%% =====================================================================
%  ROI DETECTION WITH BRAIN MASKING
%  =====================================================================
function roi_data = detect_rois_with_masking(filepath, config, logFid)

reader = bfGetReader(filepath);
reader.setSeries(0);
width = reader.getSizeX();
height = reader.getSizeY();
nFrames = reader.getImageCount();

fprintf('    Input: %dx%d, %d frames\n', height, width, nFrames);

% Create brain mask
fprintf('    Creating brain mask...\n');
[brain_mask, meanImg, stdImg] = create_brain_mask(reader, height, width, nFrames);
mask_coverage = sum(brain_mask(:)) / numel(brain_mask) * 100;
fprintf('    Brain mask: %.1f%% of FOV\n', mask_coverage);
logf(logFid, '    Brain mask: %.1f%%\n', mask_coverage);

% Calculate Cn/PNR maps
fprintf('    Calculating Cn/PNR maps...\n');
[Cn, PNR] = calculate_maps_masked(reader, height, width, nFrames, brain_mask);

% Detect ROIs
fprintf('    Finding ROI seeds (v4.4 corrected sizes)...\n');
[roi_masks, roi_centers, nROIs, detection_info] = detect_rois_v44(Cn, PNR, meanImg, stdImg, brain_mask, config, logFid);
fprintf('    Detected %d ROIs\n', nROIs);

% Store results
roi_data = struct();
roi_data.Cn = Cn;
roi_data.PNR = PNR;
roi_data.meanImg = meanImg;
roi_data.stdImg = stdImg;
roi_data.brain_mask = brain_mask;
roi_data.dims = [height, width, nFrames];
roi_data.detection_info = detection_info;

if nROIs == 0
    fprintf('    ⚠ WARNING: No ROIs detected!\n');
    logf(logFid, '    ⚠ No ROIs detected\n');
    roi_data.roi_masks = [];
    roi_data.roi_centers = [];
    roi_data.F = [];
    roi_data.Fneu = [];
    roi_data.dFF = [];
    roi_data.snr = [];
    roi_data.event_rate = [];
    roi_data.event_amplitude = [];
    roi_data.event_count = [];
    reader.close();
    return;
end

% Extract traces
fprintf('    Extracting traces...\n');
[F, Fneu, dFF, dFF_raw, quality_metrics] = extract_traces_full(reader, roi_masks, brain_mask, config);

% Apply Savitzky-Golay smoothing
fprintf('    Applying Savitzky-Golay smoothing...\n');
dFF_smooth = apply_savitzky_golay(dFF, config.frameRate);

% QC
fprintf('    Applying QC...\n');
[keep_idx, snr, event_rate, event_amplitude, event_count] = compute_qc_metrics(dFF_smooth, config);
fprintf('    %d/%d ROIs pass QC (SNR >= %.1f)\n', sum(keep_idx), nROIs, config.roi.min_snr);
logf(logFid, '    QC: %d/%d pass\n', sum(keep_idx), nROIs);

% Filter ROIs
roi_data.roi_masks = roi_masks(:,:,keep_idx);
roi_data.roi_centers = roi_centers(keep_idx,:);
roi_data.F = F(keep_idx,:);
roi_data.Fneu = Fneu(keep_idx,:);
roi_data.dFF = dFF(keep_idx,:);
roi_data.dFF_raw = dFF_raw(keep_idx,:);  % v4.5: Before global regression
roi_data.dFF_smooth = dFF_smooth(keep_idx,:);
roi_data.snr = snr(keep_idx);
roi_data.event_rate = event_rate(keep_idx);
roi_data.event_amplitude = event_amplitude(keep_idx);
roi_data.event_count = event_count(keep_idx);

% v4.5: Store quality metrics
roi_data.quality_metrics = quality_metrics;
if isfield(quality_metrics, 'global_corr_per_roi_before')
    roi_data.quality_metrics.global_corr_per_roi_before = quality_metrics.global_corr_per_roi_before(keep_idx);
end
if isfield(quality_metrics, 'global_corr_per_roi_after')
    roi_data.quality_metrics.global_corr_per_roi_after = quality_metrics.global_corr_per_roi_after(keep_idx);
end

% Synchrony detection
if sum(keep_idx) >= 2
    [sync_score, sync_matrix] = compute_synchrony_score(dFF_smooth(keep_idx,:));
    roi_data.synchrony_score = sync_score;
    roi_data.synchrony_matrix = sync_matrix;
    
    if sync_score > 0.5
        fprintf('    ⚠ HIGH SYNCHRONY (%.2f) - possible motion artifact!\n', sync_score);
    elseif sync_score > 0.3
        fprintf('    → Moderate synchrony (%.2f)\n', sync_score);
    else
        fprintf('    ✓ Low synchrony (%.2f) - traces appear independent\n', sync_score);
    end
else
    roi_data.synchrony_score = NaN;
    roi_data.synchrony_matrix = [];
end

if isfield(detection_info, 'roi_levels')
    roi_data.detection_info.roi_levels = detection_info.roi_levels(keep_idx);
end

reader.close();

end

%% =====================================================================
%  SAVITZKY-GOLAY SMOOTHING
%  =====================================================================
function dFF_smooth = apply_savitzky_golay(dFF, frameRate)

nROIs = size(dFF, 1);
dFF_smooth = zeros(size(dFF), 'single');

window_frames = round(0.15 * frameRate);
if mod(window_frames, 2) == 0
    window_frames = window_frames + 1;
end
window_frames = max(5, window_frames);

poly_order = min(3, window_frames - 1);

for i = 1:nROIs
    trace = dFF(i, :);
    try
        dFF_smooth(i, :) = sgolayfilt(trace, poly_order, window_frames);
    catch
        dFF_smooth(i, :) = smoothdata(trace, 'movmean', window_frames);
    end
end

end

%% =====================================================================
%  SYNCHRONY DETECTION
%  =====================================================================
function [sync_score, sync_matrix] = compute_synchrony_score(dFF)

nROIs = size(dFF, 1);

if nROIs < 2
    sync_score = NaN;
    sync_matrix = [];
    return;
end

sync_matrix = corrcoef(dFF');
upper_tri = triu(true(nROIs), 1);
pairwise_corrs = sync_matrix(upper_tri);
sync_score = mean(abs(pairwise_corrs), 'omitnan');

end

%% =====================================================================
%  BRAIN MASK CREATION
%  =====================================================================
function [brain_mask, meanImg, stdImg] = create_brain_mask(reader, height, width, nFrames)

nSample = min(200, nFrames);
sample_idx = round(linspace(1, nFrames, nSample));

meanImg = zeros(height, width, 'single');
M2 = zeros(height, width, 'single');

for i = 1:length(sample_idx)
    img = single(bfGetPlane(reader, sample_idx(i)));
    delta = img - meanImg;
    meanImg = meanImg + delta / i;
    M2 = M2 + delta .* (img - meanImg);
end
stdImg = sqrt(M2 / (length(sample_idx) - 1));

mean_norm = (meanImg - min(meanImg(:))) / (max(meanImg(:)) - min(meanImg(:)) + eps);
std_norm = (stdImg - min(stdImg(:))) / (max(stdImg(:)) - min(stdImg(:)) + eps);

combined = imgaussfilt(mean_norm .* std_norm, 10);
threshold = graythresh(combined) * 0.5;
brain_mask = combined > threshold;

brain_mask = imfill(brain_mask, 'holes');
brain_mask = bwareaopen(brain_mask, 5000);

cc = bwconncomp(brain_mask);
if cc.NumObjects > 0
    areas = cellfun(@numel, cc.PixelIdxList);
    [~, largest] = max(areas);
    brain_mask = false(size(brain_mask));
    brain_mask(cc.PixelIdxList{largest}) = true;
end

brain_mask = imdilate(brain_mask, strel('disk', 5));

end

%% =====================================================================
%  CALCULATE Cn AND PNR MAPS
%  =====================================================================
function [Cn, PNR] = calculate_maps_masked(reader, height, width, nFrames, brain_mask)

nSample = min(500, nFrames);
sample_idx = round(linspace(1, nFrames, nSample));

chunk_size = 50;
Cn = zeros(height, width, 'single');
PNR_accum = zeros(height, width, 'single');
n_chunks = 0;

for start_i = 1:chunk_size:nSample
    end_i = min(start_i + chunk_size - 1, nSample);
    chunk_frames = end_i - start_i + 1;
    
    Y = zeros(height, width, chunk_frames, 'single');
    for j = 1:chunk_frames
        img = single(bfGetPlane(reader, sample_idx(start_i + j - 1)));
        img(~brain_mask) = NaN;
        Y(:,:,j) = img;
    end
    
    Y_mean = nanmean(Y, 3);
    Y_std = nanstd(Y, 0, 3);
    Y_norm = (Y - Y_mean) ./ (Y_std + eps);
    Y_norm(isnan(Y_norm)) = 0;
    
    Cn_chunk = zeros(height, width);
    for dy = -1:1
        for dx = -1:1
            if dx == 0 && dy == 0, continue; end
            Y_shift = circshift(Y_norm, [dy, dx, 0]);
            Cn_chunk = Cn_chunk + mean(Y_norm .* Y_shift, 3);
        end
    end
    Cn = Cn + Cn_chunk / 8;
    
    Y_max = max(Y, [], 3);
    PNR_accum = PNR_accum + (Y_max - Y_mean) ./ (Y_std + eps);
    
    n_chunks = n_chunks + 1;
end

Cn = Cn / n_chunks;
PNR = PNR_accum / n_chunks;

Cn(~brain_mask) = 0;
PNR(~brain_mask) = 0;

Cn = imgaussfilt(Cn, 1.5);
PNR = imgaussfilt(PNR, 1);

end

%% =====================================================================
%  ROI DETECTION v4.4
%  =====================================================================
function [roi_masks, roi_centers, nROIs, detection_info] = detect_rois_v44(Cn, PNR, meanImg, stdImg, brain_mask, config, logFid)

[h, w] = size(Cn);
detection_info = struct();

Cn_brain = Cn .* brain_mask;
PNR_brain = PNR .* brain_mask;

cn_vals = Cn_brain(brain_mask);
pnr_vals = PNR_brain(brain_mask);

% Store diagnostics
detection_info.cn_mean = mean(cn_vals);
detection_info.cn_std = std(cn_vals);
detection_info.cn_median = median(cn_vals);
detection_info.cn_max = max(cn_vals);
detection_info.cn_min = min(cn_vals);
detection_info.cn_range = max(cn_vals) - min(cn_vals);
detection_info.cn_prctiles = prctile(cn_vals, [25, 50, 75, 90, 95]);
detection_info.pnr_mean = mean(pnr_vals);
detection_info.pnr_std = std(pnr_vals);
detection_info.pnr_median = median(pnr_vals);
detection_info.pnr_max = max(pnr_vals);
detection_info.pnr_min = min(pnr_vals);
detection_info.pnr_prctiles = prctile(pnr_vals, [25, 50, 75, 90, 95]);

fprintf('    Cn: median=%.3f, range=%.3f, max=%.3f\n', ...
    detection_info.cn_median, detection_info.cn_range, detection_info.cn_max);
fprintf('    PNR: median=%.2f, 75th=%.2f, max=%.2f\n', ...
    detection_info.pnr_median, detection_info.pnr_prctiles(3), detection_info.pnr_max);

% Check if Cn is uniformly high
cn_is_uniform = detection_info.cn_range < 0.15 && detection_info.cn_median > 0.7;
if cn_is_uniform
    fprintf('    ⚠ Cn is uniformly high - using PNR-only detection\n');
    detection_info.cn_uniform = true;
else
    detection_info.cn_uniform = false;
end

% Check data quality
if detection_info.pnr_max < config.roi.min_pnr
    detection_info.data_quality = 'low_pnr';
elseif detection_info.cn_max < config.roi.min_corr
    detection_info.data_quality = 'poor_correlation';
else
    detection_info.data_quality = 'normal';
end

% Weak signal mode
weak_signal_mode = false;
effective_min_pnr = config.roi.min_pnr;
effective_min_corr = config.roi.min_corr;

if detection_info.pnr_prctiles(3) < config.roi.min_pnr
    weak_signal_mode = true;
    effective_min_pnr = config.roi.min_pnr * 0.75;
    effective_min_corr = config.roi.min_corr * 0.80;
    fprintf('    → WEAK SIGNAL MODE: Lowering thresholds\n');
end

detection_info.weak_signal_mode = weak_signal_mode;
detection_info.effective_min_pnr = effective_min_pnr;
detection_info.effective_min_corr = effective_min_corr;

% Multi-scale seed detection
scales = [config.roi.gSig * 0.5, config.roi.gSig, config.roi.gSig * 1.5];
all_seeds_y = [];
all_seeds_x = [];
all_seeds_score = [];
all_seeds_level = [];

fprintf('    Searching for seeds...\n');

% LEVEL 1-3: Cn/PNR detection
if ~cn_is_uniform
    cn_adaptive = prctile(cn_vals, 75);
    pnr_adaptive = prctile(pnr_vals, 75);
    
    cn_thresh_strict = max(effective_min_corr, min(cn_adaptive, 0.6));
    pnr_thresh_strict = max(effective_min_pnr, min(pnr_adaptive, 6.0));
    
    threshold_levels = {
        struct('cn', cn_thresh_strict, 'pnr', pnr_thresh_strict, 'name', 'strict');
        struct('cn', max(effective_min_corr, cn_thresh_strict * 0.7), ...
               'pnr', max(effective_min_pnr * 0.8, pnr_thresh_strict * 0.7), 'name', 'medium');
        struct('cn', effective_min_corr, 'pnr', effective_min_pnr, 'name', 'relaxed');
    };
    
    detection_info.threshold_levels = threshold_levels;
    detection_info.seeds_per_level = zeros(3, 1);
    
    for level_idx = 1:length(threshold_levels)
        level = threshold_levels{level_idx};
        seeds_this_level = 0;
        
        seed_mask_base = (Cn_brain > level.cn) & (PNR_brain > level.pnr);
        
        for scale = scales
            Cn_smooth = imgaussfilt(Cn_brain, scale);
            local_max = imregionalmax(Cn_smooth);
            seed_mask = seed_mask_base & local_max & brain_mask;
            
            [y_new, x_new] = find(seed_mask);
            
            if ~isempty(y_new)
                scores_new = arrayfun(@(yi,xi) Cn(yi,xi) * PNR(yi,xi), y_new, x_new);
                
                all_seeds_y = [all_seeds_y; y_new];
                all_seeds_x = [all_seeds_x; x_new];
                all_seeds_score = [all_seeds_score; scores_new];
                all_seeds_level = [all_seeds_level; repmat(level_idx, length(y_new), 1)];
                seeds_this_level = seeds_this_level + length(y_new);
            end
        end
        
        detection_info.seeds_per_level(level_idx) = seeds_this_level;
        fprintf('      %s (Cn>%.2f, PNR>%.1f): %d seeds\n', level.name, level.cn, level.pnr, seeds_this_level);
    end
else
    detection_info.threshold_levels = {};
    detection_info.seeds_per_level = [];
end

% LEVEL 4: PNR-only detection
if isempty(all_seeds_y) || cn_is_uniform
    fprintf('    Trying PNR-only detection...\n');
    
    pnr_thresh = max(effective_min_pnr * 0.9, prctile(pnr_vals, 70));
    
    for scale = scales
        PNR_smooth = imgaussfilt(PNR_brain, scale);
        local_max = imregionalmax(PNR_smooth);
        pnr_seeds = local_max & (PNR_brain > pnr_thresh) & brain_mask;
        
        [y_new, x_new] = find(pnr_seeds);
        
        if ~isempty(y_new)
            scores_new = arrayfun(@(yi,xi) PNR(yi,xi), y_new, x_new);
            
            all_seeds_y = [all_seeds_y; y_new];
            all_seeds_x = [all_seeds_x; x_new];
            all_seeds_score = [all_seeds_score; scores_new];
            all_seeds_level = [all_seeds_level; repmat(4, length(y_new), 1)];
        end
    end
    
    n_pnr_seeds = sum(all_seeds_level == 4);
    fprintf('      PNR-only: %d seeds\n', n_pnr_seeds);
    detection_info.n_seeds_pnr_only = n_pnr_seeds;
end

% LEVEL 5-6: Fallback methods
if isempty(all_seeds_y)
    fprintf('    Trying fallback methods...\n');
    
    mean_brain = meanImg .* brain_mask;
    mean_smooth = imgaussfilt(mean_brain, config.roi.gSig);
    local_max = imregionalmax(mean_smooth);
    intensity_thresh = prctile(mean_smooth(brain_mask), 40);
    fallback_seeds = local_max & (mean_smooth > intensity_thresh) & brain_mask;
    
    [y_new, x_new] = find(fallback_seeds);
    
    if ~isempty(y_new)
        scores_new = arrayfun(@(yi,xi) mean_smooth(yi,xi), y_new, x_new);
        all_seeds_y = [all_seeds_y; y_new];
        all_seeds_x = [all_seeds_x; x_new];
        all_seeds_score = [all_seeds_score; scores_new];
        all_seeds_level = [all_seeds_level; repmat(5, length(y_new), 1)];
    end
end

% Remove duplicates and apply NMS
if ~isempty(all_seeds_y)
    [~, unique_idx] = unique([all_seeds_y, all_seeds_x], 'rows', 'first');
    all_seeds_y = all_seeds_y(unique_idx);
    all_seeds_x = all_seeds_x(unique_idx);
    all_seeds_score = all_seeds_score(unique_idx);
    all_seeds_level = all_seeds_level(unique_idx);
    
    [~, sort_idx] = sort(all_seeds_score, 'descend');
    all_seeds_y = all_seeds_y(sort_idx);
    all_seeds_x = all_seeds_x(sort_idx);
    all_seeds_score = all_seeds_score(sort_idx);
    all_seeds_level = all_seeds_level(sort_idx);
    
    % NMS with v4.4 reduced spacing
    min_seed_spacing = config.roi.gSig * 1.5;
    keep_seed = true(length(all_seeds_y), 1);
    
    for i = 1:length(all_seeds_y)
        if ~keep_seed(i), continue; end
        for j = (i+1):length(all_seeds_y)
            if ~keep_seed(j), continue; end
            dist = sqrt((all_seeds_y(i) - all_seeds_y(j))^2 + (all_seeds_x(i) - all_seeds_x(j))^2);
            if dist < min_seed_spacing
                keep_seed(j) = false;
            end
        end
    end
    
    all_seeds_y = all_seeds_y(keep_seed);
    all_seeds_x = all_seeds_x(keep_seed);
    all_seeds_score = all_seeds_score(keep_seed);
    all_seeds_level = all_seeds_level(keep_seed);
    
    fprintf('    After NMS: %d seeds\n', length(all_seeds_y));
end

detection_info.n_seeds = length(all_seeds_y);

if isempty(all_seeds_y)
    fprintf('    ⚠ All detection methods failed\n');
    roi_masks = [];
    roi_centers = [];
    nROIs = 0;
    detection_info.n_rois_detected = 0;
    detection_info.roi_levels = [];
    detection_info.n_rois_by_level = zeros(1, 6);
    detection_info.failure_reason = 'no_seeds_found';
    return;
end

% Grow ROIs with v4.4 corrected sizes
maxROIs = min(length(all_seeds_y), config.roi.max_rois);
roi_masks = false(h, w, maxROIs);
roi_centers = zeros(maxROIs, 2);
roi_detection_level = zeros(maxROIs, 1);

used = false(h, w);
count = 0;

target_area = pi * (config.roi.gSig * config.roi.target_area_factor)^2;
target_area = min(target_area, config.roi.max_area_px);

fprintf('    Growing ROIs (target_area=%.0f px, max=%.0f px)...\n', target_area, config.roi.max_area_px);

for i = 1:length(all_seeds_y)
    if count >= maxROIs, break; end
    
    y = all_seeds_y(i);
    x = all_seeds_x(i);
    seed_level = all_seeds_level(i);
    
    if used(y, x), continue; end
    
    roi_mask = grow_roi_local_adaptive_v44(Cn_brain, PNR_brain, meanImg, stdImg, ...
        x, y, config, brain_mask, used, h, w, seed_level, target_area);
    
    area = sum(roi_mask(:));
    
    if area < config.roi.min_area_px || area > config.roi.max_area_px
        continue;
    end
    
    if any(used(roi_mask))
        overlap_frac = sum(used(roi_mask)) / area;
        if overlap_frac > 0.3, continue; end
        roi_mask = roi_mask & ~used;
        area = sum(roi_mask(:));
        if area < config.roi.min_area_px, continue; end
    end
    
    count = count + 1;
    roi_masks(:,:,count) = roi_mask;
    roi_centers(count,:) = [x, y];
    roi_detection_level(count) = seed_level;
    used(roi_mask) = true;
end

roi_masks = roi_masks(:,:,1:count);
roi_centers = roi_centers(1:count,:);
roi_detection_level = roi_detection_level(1:count);
nROIs = count;

detection_info.n_rois_detected = nROIs;
detection_info.roi_levels = roi_detection_level;
detection_info.n_rois_by_level = zeros(1, 6);
for lev = 1:6
    detection_info.n_rois_by_level(lev) = sum(roi_detection_level == lev);
end

fprintf('    Final: %d ROIs\n', nROIs);

end

%% =====================================================================
%  v4.4: ROI GROWING (with corrected size constraints)
%  =====================================================================
function roi_mask = grow_roi_local_adaptive_v44(Cn_brain, PNR_brain, meanImg, stdImg, x, y, config, brain_mask, used, h, w, seed_level, target_area)

gSig = config.roi.gSig;
[xx, yy] = meshgrid(1:w, 1:h);

% Robust seed values
y_lo = max(1, y-1); y_hi = min(h, y+1);
x_lo = max(1, x-1); x_hi = min(w, x+1);

seed_region_cn = Cn_brain(y_lo:y_hi, x_lo:x_hi);
seed_region_pnr = PNR_brain(y_lo:y_hi, x_lo:x_hi);
seed_region_mean = meanImg(y_lo:y_hi, x_lo:x_hi);
seed_region_std = stdImg(y_lo:y_hi, x_lo:x_hi);

seed_cn = max(seed_region_cn(:));
seed_pnr = max(seed_region_pnr(:));
seed_mean = max(seed_region_mean(:));
seed_std = max(seed_region_std(:));

% v4.4: Tighter search region for smaller neurons
search_half = round(gSig * 2.0);  % Reduced from 2.5
y_min = max(1, y - search_half);
y_max = min(h, y + search_half);
x_min = max(1, x - search_half);
x_max = min(w, x + search_half);

search_region = false(h, w);
search_region(y_min:y_max, x_min:x_max) = true;
search_region = search_region & brain_mask & ~used;

if sum(search_region(:)) < 5
    distance_map = sqrt((xx - x).^2 + (yy - y).^2);
    roi_mask = (distance_map <= gSig) & brain_mask & ~used;
    return;
end

% Seed-relative thresholds
switch seed_level
    case {1, 2, 3}
        cn_thresh = seed_cn * 0.80;
        pnr_thresh = seed_pnr * 0.70;
        grow_mask = (Cn_brain >= cn_thresh) & (PNR_brain >= pnr_thresh) & search_region;
    case 4
        pnr_thresh = seed_pnr * 0.75;
        grow_mask = (PNR_brain >= pnr_thresh) & search_region;
    case 5
        mean_thresh = seed_mean * 0.85;
        grow_mask = (meanImg >= mean_thresh) & search_region;
    case 6
        std_thresh = seed_std * 0.80;
        grow_mask = (stdImg >= std_thresh) & search_region;
    otherwise
        cn_thresh = seed_cn * 0.80;
        grow_mask = (Cn_brain >= cn_thresh) & search_region;
end

roi_mask = get_seed_component(grow_mask, x, y, h, w);
area = sum(roi_mask(:));

% Relaxed thresholds if too small
if area < config.roi.min_area_px
    switch seed_level
        case {1, 2, 3}
            cn_thresh = seed_cn * 0.65;
            pnr_thresh = seed_pnr * 0.55;
            grow_mask = (Cn_brain >= cn_thresh) & (PNR_brain >= pnr_thresh) & search_region;
        case 4
            pnr_thresh = seed_pnr * 0.60;
            grow_mask = (PNR_brain >= pnr_thresh) & search_region;
        otherwise
            distance_map = sqrt((xx - x).^2 + (yy - y).^2);
            grow_mask = (distance_map <= gSig * 1.2) & brain_mask & ~used;
    end
    
    roi_mask = get_seed_component(grow_mask, x, y, h, w);
    area = sum(roi_mask(:));
end

% Circular fallback
if area < config.roi.min_area_px
    distance_map = sqrt((xx - x).^2 + (yy - y).^2);
    roi_mask = (distance_map <= gSig * 1.2) & brain_mask & ~used;
    area = sum(roi_mask(:));
end

% Handle oversized ROIs with v4.4 stricter limits
if area > config.roi.max_area_px
    % Tighten threshold
    switch seed_level
        case {1, 2, 3}
            cn_thresh = seed_cn * 0.92;
            pnr_thresh = seed_pnr * 0.88;
            grow_mask = (Cn_brain >= cn_thresh) & (PNR_brain >= pnr_thresh) & search_region;
        case 4
            pnr_thresh = seed_pnr * 0.92;
            grow_mask = (PNR_brain >= pnr_thresh) & search_region;
        otherwise
            distance_map = sqrt((xx - x).^2 + (yy - y).^2);
            grow_mask = (distance_map <= gSig) & brain_mask & ~used;
    end
    
    roi_mask_tight = get_seed_component(grow_mask, x, y, h, w);
    area_tight = sum(roi_mask_tight(:));
    
    if area_tight >= config.roi.min_area_px && area_tight <= config.roi.max_area_px
        roi_mask = roi_mask_tight;
    else
        % Morphological erosion
        se = strel('disk', 1);
        roi_mask_temp = roi_mask;
        while sum(roi_mask_temp(:)) > config.roi.max_area_px
            roi_mask_eroded = imerode(roi_mask_temp, se);
            if roi_mask_eroded(y, x)
                roi_mask_temp = roi_mask_eroded;
            else
                break;
            end
        end
        if sum(roi_mask_temp(:)) >= config.roi.min_area_px
            roi_mask = roi_mask_temp;
        else
            distance_map = sqrt((xx - x).^2 + (yy - y).^2);
            target_radius = sqrt(target_area / pi);
            roi_mask = (distance_map <= target_radius) & brain_mask & ~used;
        end
    end
end

% Ensure ROI contains seed
if ~roi_mask(y, x)
    distance_map = sqrt((xx - x).^2 + (yy - y).^2);
    roi_mask = (distance_map <= gSig) & brain_mask & ~used;
end

end

%% =====================================================================
%  HELPER: Get connected component containing seed
%  =====================================================================
function roi_mask = get_seed_component(grow_mask, x, y, h, w)
    CC = bwconncomp(grow_mask, 8);
    label_map = labelmatrix(CC);
    
    if y >= 1 && y <= h && x >= 1 && x <= w
        seed_label = label_map(y, x);
    else
        seed_label = 0;
    end
    
    if seed_label > 0
        roi_mask = false(h, w);
        roi_mask(CC.PixelIdxList{seed_label}) = true;
    else
        roi_mask = false(h, w);
    end
end

%% =====================================================================
%  TRACE EXTRACTION
%  =====================================================================
function [F, Fneu, dFF, dFF_raw, quality_metrics] = extract_traces_full(reader, roi_masks, brain_mask, config)
% EXTRACT_TRACES_FULL - Extract fluorescence traces with optional global regression
%
% v4.5 UPDATE: Added global signal regression for widefield imaging
%   - dFF_raw: dF/F before global regression (standard neuropil correction only)
%   - dFF: dF/F after global regression (cleaned signal)
%   - quality_metrics: struct with synchrony, global_correlation, etc.
%
% INPUTS:
%   reader     - BioFormats reader object
%   roi_masks  - [H x W x nROIs] logical masks
%   brain_mask - [H x W] logical mask of brain region
%   config     - Configuration struct

nROIs = size(roi_masks, 3);
nFrames = reader.getImageCount();

F = zeros(nROIs, nFrames, 'single');
Fneu = zeros(nROIs, nFrames, 'single');

neuropil_masks = false(size(roi_masks));
for i = 1:nROIs
    % v4.4: Ensure minimum neuropil annulus width of 4 pixels
    inner_radius = max(round(config.roi.gSig), 2);
    outer_radius = max(round(config.roi.gSig * 2.5), inner_radius + 4);  % At least 4px wide
    
    dilated = imdilate(roi_masks(:,:,i), strel('disk', outer_radius));
    inner = imdilate(roi_masks(:,:,i), strel('disk', inner_radius));
    neuropil_masks(:,:,i) = dilated & ~inner;
    
    for j = 1:nROIs
        if j ~= i
            neuropil_masks(:,:,i) = neuropil_masks(:,:,i) & ~roi_masks(:,:,j);
        end
    end
end

fprintf('    ');
for frame = 1:nFrames
    if mod(frame, 1000) == 0, fprintf('.'); end
    
    img = single(bfGetPlane(reader, frame));
    
    for roi = 1:nROIs
        F(roi, frame) = mean(img(roi_masks(:,:,roi)));
        neu_mask = neuropil_masks(:,:,roi);
        if any(neu_mask(:))
            Fneu(roi, frame) = median(img(neu_mask));
        end
    end
end
fprintf('\n');

% Calculate dF/F with standard neuropil subtraction
dFF_raw = zeros(size(F), 'single');
for roi = 1:nROIs
    F_corrected = F(roi,:) - config.roi.neuropil_coef * Fneu(roi,:);
    
    win_size = round(30 * config.frameRate);
    baseline = movmin(F_corrected, win_size);
    baseline = movmean(baseline, win_size);
    
    if mean(baseline) <= 0
        baseline = prctile(F_corrected, 8);
    end
    
    if mean(baseline) > 0
        % v4.4: Add eps to prevent division by zero in edge cases
        dFF_raw(roi,:) = (F_corrected - baseline) ./ (baseline + eps);
    end
end

% ═══════════════════════════════════════════════════════════════════════════
% v4.5 NEW: GLOBAL SIGNAL REGRESSION FOR WIDEFIELD IMAGING
% ═══════════════════════════════════════════════════════════════════════════
% This removes shared neuropil/background signal that affects all ROIs
% Critical for widefield microscopy which lacks optical sectioning
%
% IMPORTANT NOTES (per Opus 4.6 review):
% 1. Global signal is computed from brain mask EXCLUDING ROIs (true neuropil)
% 2. Linear regression includes intercept (consistent with robust method)
% 3. Validation uses pairwise synchrony, not circular post-regression corr
% 4. Both dFF_raw and dFF are saved so you can compare with/without regression

quality_metrics = struct();
quality_metrics.global_regression_applied = false;

% ─────────────────────────────────────────────────────────────────────────
% COMPUTE SYNCHRONY (PAIRWISE CORRELATION) BEFORE REGRESSION
% This is the proper metric - not correlation with mean (which is circular)
% ─────────────────────────────────────────────────────────────────────────
if nROIs >= 2
    corr_matrix_before = corrcoef(dFF_raw');
    upper_tri = triu(true(nROIs), 1);
    pairwise_corrs_before = corr_matrix_before(upper_tri);
    quality_metrics.synchrony_before = mean(abs(pairwise_corrs_before), 'omitnan');
    quality_metrics.pairwise_corr_matrix_before = corr_matrix_before;
else
    quality_metrics.synchrony_before = NaN;
    quality_metrics.pairwise_corr_matrix_before = [];
end

% Also compute correlation with ROI-mean for backwards compatibility
global_signal_roi_mean = mean(dFF_raw, 1);
corrs_with_mean = zeros(nROIs, 1);
for i = 1:nROIs
    cc = corrcoef(dFF_raw(i,:), global_signal_roi_mean);
    if numel(cc) >= 2
        corrs_with_mean(i) = cc(1,2);
    end
end
quality_metrics.global_corr_before = mean(corrs_with_mean, 'omitnan');
quality_metrics.global_corr_per_roi_before = corrs_with_mean;

% Check for artifacts (extreme dF/F values)
max_dff = max(dFF_raw(:));
min_dff = min(dFF_raw(:));
quality_metrics.max_dff = max_dff;
quality_metrics.min_dff = min_dff;

if max_dff > 5 || min_dff < -2
    quality_metrics.artifact_flag = true;
    quality_metrics.artifact_reason = sprintf('Extreme dF/F: max=%.1f, min=%.1f', max_dff, min_dff);
    fprintf('    ⚠ ARTIFACT DETECTED: %s\n', quality_metrics.artifact_reason);
else
    quality_metrics.artifact_flag = false;
    quality_metrics.artifact_reason = '';
end

% ─────────────────────────────────────────────────────────────────────────
% APPLY GLOBAL SIGNAL REGRESSION IF ENABLED
% ─────────────────────────────────────────────────────────────────────────
if isfield(config.roi, 'apply_global_regression') && config.roi.apply_global_regression && ~quality_metrics.artifact_flag
    fprintf('    Applying global signal regression...\n');
    
    % ─────────────────────────────────────────────────────────────────────
    % COMPUTE TRUE NEUROPIL SIGNAL FROM BRAIN MASK EXCLUDING ROIs
    % This is more accurate than using ROI mean (which is spatially biased)
    % ─────────────────────────────────────────────────────────────────────
    
    % We need to read frames again to get the true neuropil signal
    % For efficiency, sample a subset of frames
    n_sample_frames = min(500, nFrames);
    sample_idx = round(linspace(1, nFrames, n_sample_frames));
    
    % Create mask of brain excluding all ROIs
    % brain_mask is passed as parameter to this function
    roi_union = false(size(brain_mask));
    for i = 1:nROIs
        roi_union = roi_union | roi_masks(:,:,i);
    end
    neuropil_only_mask = brain_mask & ~roi_union;
    n_neuropil_pixels = sum(neuropil_only_mask(:));
    
    if n_neuropil_pixels > 100
        % Compute true neuropil signal from non-ROI brain pixels
        fprintf('    Computing neuropil reference from %d non-ROI pixels...\n', n_neuropil_pixels);
        
        % For full frames, we'll approximate: scale the sampled neuropil signal
        % to match full trace timing via interpolation
        neuropil_sampled = zeros(1, n_sample_frames);
        for idx = 1:n_sample_frames
            frame = sample_idx(idx);
            img = single(bfGetPlane(reader, frame));
            neuropil_sampled(idx) = mean(img(neuropil_only_mask));
        end
        
        % Interpolate to full frame count
        global_signal = interp1(sample_idx, neuropil_sampled, 1:nFrames, 'pchip');
        
        % Normalize to match dF/F scale (subtract baseline, divide by baseline)
        win_size = round(30 * config.frameRate);
        neu_baseline = movmin(global_signal, min(win_size, length(global_signal)));
        neu_baseline = movmean(neu_baseline, min(win_size, length(global_signal)));
        if mean(neu_baseline) > 0
            global_signal = (global_signal - neu_baseline) ./ (neu_baseline + eps);
        end
        
        quality_metrics.neuropil_source = 'brain_mask_excluding_rois';
    else
        % Fallback to ROI mean if not enough neuropil pixels
        fprintf('    Warning: Only %d neuropil pixels, using ROI mean as fallback\n', n_neuropil_pixels);
        global_signal = global_signal_roi_mean;
        quality_metrics.neuropil_source = 'roi_mean_fallback';
    end
    
    % ─────────────────────────────────────────────────────────────────────
    % REGRESS OUT GLOBAL SIGNAL FROM EACH ROI
    % Now with proper intercept in both standard and robust methods
    % ─────────────────────────────────────────────────────────────────────
    dFF = zeros(size(dFF_raw), 'single');
    
    for roi = 1:nROIs
        trace = dFF_raw(roi, :);
        
        if strcmpi(config.roi.global_regression_method, 'robust')
            % Robust regression with intercept (less sensitive to outliers)
            try
                beta = robustfit(global_signal', trace', 'bisquare');
                dFF(roi, :) = trace - (beta(1) + beta(2) * global_signal);
            catch
                % Fallback to standard with intercept
                [beta, intercept] = regress_with_intercept(global_signal, trace);
                dFF(roi, :) = trace - (intercept + beta * global_signal);
            end
        else
            % Standard linear regression WITH INTERCEPT
            % trace = intercept + beta * global + residual
            [beta, intercept] = regress_with_intercept(global_signal, trace);
            dFF(roi, :) = trace - (intercept + beta * global_signal);
        end
    end
    
    % ─────────────────────────────────────────────────────────────────────
    % COMPUTE SYNCHRONY AFTER REGRESSION (proper validation)
    % ─────────────────────────────────────────────────────────────────────
    if nROIs >= 2
        corr_matrix_after = corrcoef(dFF');
        pairwise_corrs_after = corr_matrix_after(upper_tri);
        quality_metrics.synchrony_after = mean(abs(pairwise_corrs_after), 'omitnan');
        quality_metrics.pairwise_corr_matrix_after = corr_matrix_after;
    else
        quality_metrics.synchrony_after = NaN;
        quality_metrics.pairwise_corr_matrix_after = [];
    end
    
    % Also compute ROI-mean correlation for backwards compatibility
    global_signal_after = mean(dFF, 1);
    corrs_after = zeros(nROIs, 1);
    for i = 1:nROIs
        cc = corrcoef(dFF(i,:), global_signal_after);
        if numel(cc) >= 2
            corrs_after(i) = cc(1,2);
        end
    end
    quality_metrics.global_corr_after = mean(corrs_after, 'omitnan');
    quality_metrics.global_corr_per_roi_after = corrs_after;
    quality_metrics.global_regression_applied = true;
    
    % Report using SYNCHRONY (the proper metric)
    fprintf('    Global regression: synchrony %.3f → %.3f\n', ...
        quality_metrics.synchrony_before, quality_metrics.synchrony_after);
else
    % No regression - use raw dFF
    dFF = dFF_raw;
    quality_metrics.synchrony_after = quality_metrics.synchrony_before;
    quality_metrics.global_corr_after = quality_metrics.global_corr_before;
    quality_metrics.global_corr_per_roi_after = corrs_with_mean;
    quality_metrics.neuropil_source = 'none';
    
    if quality_metrics.artifact_flag
        fprintf('    Skipping global regression due to artifact flag\n');
    elseif ~isfield(config.roi, 'apply_global_regression') || ~config.roi.apply_global_regression
        fprintf('    Global regression disabled\n');
    end
end

% Quality assessment based on synchrony (more meaningful than ROI-mean corr)
if quality_metrics.synchrony_before > 0.5
    quality_metrics.neuropil_severity = 'SEVERE';
elseif quality_metrics.synchrony_before > 0.3
    quality_metrics.neuropil_severity = 'HIGH';
elseif quality_metrics.synchrony_before > 0.15
    quality_metrics.neuropil_severity = 'MILD';
else
    quality_metrics.neuropil_severity = 'LOW';
end

end

%% =====================================================================
%  HELPER: Linear regression with intercept
%  =====================================================================
function [beta, intercept] = regress_with_intercept(x, y)
% Linear regression: y = intercept + beta * x
% Returns beta (slope) and intercept

x = x(:);  % Column vector
y = y(:);

% Demean for numerical stability
x_mean = mean(x);
y_mean = mean(y);
x_centered = x - x_mean;
y_centered = y - y_mean;

% Compute slope
beta = (x_centered' * y_centered) / (x_centered' * x_centered + eps);

% Compute intercept
intercept = y_mean - beta * x_mean;

end

%% =====================================================================
%  QC METRICS
%  =====================================================================
function [keep_idx, snr, event_rate, event_amplitude, event_count] = compute_qc_metrics(dFF, config)

nROIs = size(dFF, 1);
nFrames = size(dFF, 2);
fs = config.frameRate;
duration_s = nFrames / fs;

min_width_frames = max(3, round(config.events.min_width_s * fs));
max_width_frames = round(config.events.max_width_s * fs);
min_distance_frames = max(3, round(config.events.min_distance_s * fs));

snr = zeros(nROIs, 1);
event_rate = zeros(nROIs, 1);
event_amplitude = zeros(nROIs, 1);
event_count = zeros(nROIs, 1);

for i = 1:nROIs
    trace = dFF(i, :)';
    
    signal = prctile(trace, 95) - prctile(trace, 25);
    trace_smooth = smoothdata(trace, 'movmean', round(2*fs));
    noise = mad(trace - trace_smooth, 1) * 1.4826;
    snr(i) = signal / (noise + eps);
    
    win_size = min(round(30 * fs), floor(nFrames / 4));
    baseline = smoothdata(movmin(trace, win_size), 'gaussian', round(win_size / 2));
    trace_filt = medfilt1(trace - baseline, 3);
    
    try
        [pks, locs] = findpeaks(trace_filt, ...
            'MinPeakProminence', config.events.min_prominence, ...
            'MinPeakWidth', min_width_frames, ...
            'MaxPeakWidth', max_width_frames, ...
            'MinPeakDistance', min_distance_frames);
    catch
        try
            [pks, locs] = findpeaks(trace_filt, ...
                'MinPeakProminence', config.events.min_prominence, ...
                'MinPeakWidth', min_width_frames, ...
                'MinPeakDistance', min_distance_frames);
        catch
            pks = []; locs = [];
        end
    end
    
    if ~isempty(pks)
        valid = true(length(pks), 1);
        for e = 1:length(pks)
            search_start = max(1, locs(e) - round(2*fs));
            [~, min_idx] = min(trace_filt(search_start:locs(e)));
            rise_start = search_start + min_idx - 1;
            rise_time = (locs(e) - rise_start) / fs;
            if rise_time > 0
                rise_rate = (pks(e) - trace_filt(rise_start)) / rise_time;
                if rise_rate < config.events.min_rise_rate
                    valid(e) = false;
                end
            end
        end
        pks = pks(valid);
    end
    
    event_count(i) = length(pks);
    event_rate(i) = length(pks) / duration_s;
    event_amplitude(i) = mean(pks(pks > 0), 'omitnan');
    if isnan(event_amplitude(i)), event_amplitude(i) = 0; end
end

% Adaptive SNR threshold
median_snr = median(snr);
effective_min_snr = config.roi.min_snr;

if median_snr < config.roi.min_snr && median_snr > config.roi.min_snr * 0.5
    effective_min_snr = max(config.roi.min_snr * 0.7, median_snr * 0.8);
    fprintf('    → Weak signals: adjusting SNR threshold %.1f→%.1f\n', config.roi.min_snr, effective_min_snr);
end

keep_idx = snr >= effective_min_snr;
has_events = event_count >= 2;
borderline_snr = snr >= effective_min_snr * 0.8 & snr < effective_min_snr;
keep_idx = keep_idx | (has_events & borderline_snr);

fprintf('    Events: %d/%d ROIs have events, mean rate=%.3f Hz, mean SNR=%.2f\n', ...
    sum(event_count > 0), nROIs, mean(event_rate), mean(snr));

end

%% =====================================================================
%  v4.4 ENHANCED FIGURE: ANALYSIS
%  =====================================================================
function generate_analysis_figures_v45(roi_data, figDir, basename, config)

nROIs = size(roi_data.roi_centers, 1);
fig = figure('Position', [50 50 1800 1000], 'Visible', 'off');

% Color scheme for detection levels
level_colors = [
    0.2 0.7 0.2;   % L1: Green
    0.9 0.6 0.1;   % L2: Orange
    0.8 0.2 0.2;   % L3: Red
    0.5 0.2 0.8;   % L4: Purple
    0.2 0.5 0.8;   % L5: Blue
    0.2 0.8 0.8    % L6: Cyan
];

if nROIs > 0 && isfield(roi_data.detection_info, 'roi_levels')
    roi_levels = roi_data.detection_info.roi_levels;
else
    roi_levels = ones(max(nROIs, 1), 1);
end

% Panel 1: Mean image with ENHANCED ROI visualization
subplot(2,4,1);

% Build composite RGB image for better visibility
meanImg_masked = roi_data.meanImg .* roi_data.brain_mask;
bg_min = prctile(meanImg_masked(roi_data.brain_mask), 1);
bg_max = prctile(meanImg_masked(roi_data.brain_mask), 99);
bg_panel = (meanImg_masked - bg_min) / (bg_max - bg_min + eps);
bg_panel = max(0, min(1, bg_panel));
bg_panel = bg_panel .^ 0.5;  % Gamma brighten
bg_panel = bg_panel * 0.55;  % Dim for overlay visibility
[ph, pw] = size(bg_panel);
panel_rgb = repmat(bg_panel, [1 1 3]);

if nROIs > 0
    % Paint ROIs by detection level color
    for i = 1:nROIs
        level = min(max(roi_levels(i), 1), 6);
        color = level_colors(level,:);
        mask = roi_data.roi_masks(:,:,i);
        
        for c = 1:3
            ch = panel_rgb(:,:,c);
            ch(mask) = ch(mask) * 0.4 + color(c) * 0.6;
            panel_rgb(:,:,c) = ch;
        end
    end
    
    % Add bright boundaries
    for i = 1:nROIs
        perim = bwperim(roi_data.roi_masks(:,:,i));
        for c = 1:3
            ch = panel_rgb(:,:,c);
            ch(perim) = 1.0;
            panel_rgb(:,:,c) = ch;
        end
    end
end

imshow(panel_rgb); axis image; hold on;
title(sprintf('%s: %d ROIs (v4.5)', basename, nROIs), 'Interpreter', 'none', 'FontSize', 9);
colorbar;

% Panel 2: Cn map
subplot(2,4,2);
imagesc(roi_data.Cn); colormap(gca, gray); axis image;
info = roi_data.detection_info;
cn_status = '';
if isfield(info, 'cn_uniform') && info.cn_uniform
    cn_status = ' [UNIFORM]';
end
title(sprintf('Cn (med=%.2f, range=%.2f)%s', info.cn_median, info.cn_range, cn_status), 'FontSize', 9);
colorbar; caxis([0, min(1, max(roi_data.Cn(:)))]);

% Panel 3: PNR map
subplot(2,4,3);
imagesc(roi_data.PNR); colormap(gca, hot); axis image;
title(sprintf('PNR (max=%.2f, thresh=%.1f)', info.pnr_max, config.roi.min_pnr), 'FontSize', 9);
colorbar;

% Panel 4: Brain mask with ROI count
subplot(2,4,4);
imagesc(roi_data.brain_mask); colormap(gca, gray); axis image;
title(sprintf('Brain Mask (%.1f%%)', sum(roi_data.brain_mask(:))/numel(roi_data.brain_mask)*100), 'FontSize', 9);

% Panel 5: ROI size distribution
subplot(2,4,5);
if nROIs > 0
    areas = zeros(nROIs, 1);
    for i = 1:nROIs
        areas(i) = sum(roi_data.roi_masks(:,:,i), 'all');
    end
    areas_um2 = areas * config.pixelSize_um^2;
    diameters_um = 2 * sqrt(areas_um2 / pi);
    
    histogram(diameters_um, 20, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    hold on;
    xline(2, 'g--', 'LineWidth', 1.5, 'Label', 'Min');
    xline(3, 'r--', 'LineWidth', 1.5, 'Label', 'Max expected');
    xlabel('Diameter (µm)'); ylabel('Count');
    title('ROI Size Distribution', 'FontSize', 9);
else
    text(0.5, 0.5, 'No ROIs', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('ROI Sizes', 'FontSize', 9);
end

% Panel 6: PNR distribution
subplot(2,4,6);
pnr_vals = roi_data.PNR(roi_data.brain_mask);
histogram(pnr_vals, 50, 'FaceColor', [0.7 0.3 0.3], 'EdgeColor', 'none');
hold on; 
xline(config.roi.min_pnr, 'k--', 'LineWidth', 2);
xlabel('PNR'); ylabel('Count'); title(sprintf('PNR (thresh=%.1f)', config.roi.min_pnr), 'FontSize', 9);

% Panel 7: ROIs by detection level
subplot(2,4,7);
if nROIs > 0
    level_counts = zeros(1, 6);
    for lev = 1:6
        level_counts(lev) = sum(roi_levels == lev);
    end
    b = bar(1:6, level_counts);
    b.FaceColor = 'flat';
    for lev = 1:6
        b.CData(lev,:) = level_colors(lev,:);
    end
    set(gca, 'XTick', 1:6, 'XTickLabel', {'L1','L2','L3','PNR','Mean','Std'});
    ylabel('# ROIs'); title('ROIs by Detection Level', 'FontSize', 9);
else
    text(0.5, 0.5, 'No ROIs', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Detection Levels', 'FontSize', 9);
end

% Panel 8: Best traces
subplot(2,4,8);
if nROIs > 0 && ~isempty(roi_data.dFF)
    time = (0:size(roi_data.dFF,2)-1) / config.frameRate;
    [~, order] = sort(roi_data.snr, 'descend');
    show_idx = order(1:min(6, length(order)));
    offset = 0;
    for idx = show_idx'
        if isfield(roi_data, 'dFF_smooth')
            trace = roi_data.dFF_smooth(idx,:);
        else
            trace = roi_data.dFF(idx,:);
        end
        level = min(max(roi_levels(idx), 1), 6);
        plot(time, trace + offset, 'Color', level_colors(level,:), 'LineWidth', 0.8);
        hold on;
        offset = offset + max(prctile(trace,99) - prctile(trace,1), 0.15) * 1.5;
    end
    xlabel('Time (s)'); ylabel('dF/F');
    if isfield(roi_data, 'synchrony_score') && ~isnan(roi_data.synchrony_score)
        title(sprintf('Best traces (sync=%.2f)', roi_data.synchrony_score), 'FontSize', 9);
    else
        title('Best traces', 'FontSize', 9);
    end
    xlim([0, min(time(end)*1.1, 60)]);
else
    text(0.5, 0.5, sprintf('Cn range=%.3f\nPNR max=%.2f\nSeeds: %d', ...
        info.cn_range, info.pnr_max, info.n_seeds), ...
        'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 9);
    title('Detection Info', 'FontSize', 9);
end

% Summary annotation
if nROIs == 0
    summary_str = sprintf('v4.5: No ROIs | Cn range=%.3f | PNR max=%.2f | Seeds: %d', ...
        info.cn_range, info.pnr_max, info.n_seeds);
else
    sync_str = '';
    if isfield(roi_data, 'synchrony_score') && ~isnan(roi_data.synchrony_score)
        if roi_data.synchrony_score > 0.5
            sync_str = sprintf(' | ⚠SYNC=%.2f', roi_data.synchrony_score);
        else
            sync_str = sprintf(' | sync=%.2f', roi_data.synchrony_score);
        end
    end
    summary_str = sprintf('v4.5: Seeds=%d → ROIs=%d | gSig=%d min=%d max=%d%s', ...
        info.n_seeds, nROIs, config.roi.gSig, config.roi.min_area_px, config.roi.max_area_px, sync_str);
end
annotation('textbox', [0.02 0.01 0.96 0.03], 'String', summary_str, ...
    'EdgeColor', 'none', 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

saveas(fig, fullfile(figDir, [basename '_analysis_v45.png']));

if config.generate_pdfs
    try
        print(fig, fullfile(figDir, [basename '_analysis_v45.pdf']), '-dpdf', '-bestfit');
    catch
    end
end

close(fig);

end

%% =====================================================================
%  v4.4 NEW FIGURE: CLEAN ROI MAP
%  =====================================================================
function generate_roi_map_v45(roi_data, figDir, basename, config)
% Generate publication-ready ROI map focused on CLEAR OUTLINES
% Design priorities: visible boundaries, high contrast, no clutter

nROIs = size(roi_data.roi_centers, 1);
if nROIs == 0
    return;
end

[h, w] = size(roi_data.meanImg);

% ─── STEP 1: Build high-contrast background ───
% Use percentile stretch + gamma correction so brain tissue is bright & clear
meanImg_masked = single(roi_data.meanImg) .* single(roi_data.brain_mask);
brain_vals = meanImg_masked(roi_data.brain_mask);
bg_lo = prctile(brain_vals, 0.5);
bg_hi = prctile(brain_vals, 99.5);
bg = (meanImg_masked - bg_lo) / (bg_hi - bg_lo + eps);
bg = max(0, min(1, bg));
bg = bg .^ 0.4;  % Strong gamma to brighten dark tissue

% ─── STEP 2: Generate distinct bright colors ───
% Use 32-color palette with large hue steps so neighbors look different
n_palette = 32;
base_hues = mod((0:n_palette-1) * 0.618033988749895, 1);  % Golden ratio spacing
base_colors = zeros(n_palette, 3);
for k = 1:n_palette
    base_colors(k,:) = hsv2rgb([base_hues(k), 0.9, 1.0]);  % High saturation, full brightness
end

% ─── STEP 3: Build RGB composite with colored outlines ───
% Start from grayscale background
bg_rgb = repmat(bg, [1 1 3]);

% Draw every ROI boundary as a thick, bright colored outline
for i = 1:nROIs
    mask = roi_data.roi_masks(:,:,i);
    perim = bwperim(mask);
    % Dilate to 2px thick for visibility
    perim_thick = imdilate(perim, strel('disk', 1));
    
    cidx = mod(i-1, n_palette) + 1;
    color = base_colors(cidx, :);
    
    for c = 1:3
        ch = bg_rgb(:,:,c);
        ch(perim_thick) = color(c);
        bg_rgb(:,:,c) = ch;
    end
end

% ─── STEP 4: Render figure ───
fig = figure('Position', [50 50 1400 1100], 'Visible', 'off', 'Color', 'k');

imshow(bg_rgb, 'Border', 'tight');
hold on;

% Scale bar (10 µm) - bottom left
scale_bar_um = 10;
scale_bar_px = scale_bar_um / config.pixelSize_um;
sb_x = w * 0.05;
sb_y = h * 0.94;
line([sb_x, sb_x + scale_bar_px], [sb_y, sb_y], 'Color', 'w', 'LineWidth', 4);
text(sb_x + scale_bar_px/2, sb_y - h*0.02, sprintf('%d µm', scale_bar_um), ...
    'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% ROI count - bottom right
text(w * 0.95, h * 0.94, sprintf('n = %d ROIs', nROIs), ...
    'Color', 'white', 'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

title(sprintf('%s: %d ROIs (v4.5)', basename, nROIs), ...
    'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');

% Save
exportgraphics_safe(fig, fullfile(figDir, [basename '_roi_map_v45.png']), 300);
if config.generate_pdfs
    try
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fullfile(figDir, [basename '_roi_map_v45.pdf']), '-dpdf', '-bestfit');
    catch
    end
end
close(fig);

% ─── FIGURE 2: Outlines on WHITE background (for print/publication) ───
fig2 = figure('Position', [50 50 1400 1100], 'Visible', 'off', 'Color', 'w');

% White background with light gray brain mask
white_bg = ones(h, w, 3, 'single');
% Tint brain region light gray so you can see the tissue shape
for c = 1:3
    ch = white_bg(:,:,c);
    ch(roi_data.brain_mask) = 0.90;  % Light gray inside brain
    white_bg(:,:,c) = ch;
end

% Draw outlines on white
for i = 1:nROIs
    mask = roi_data.roi_masks(:,:,i);
    perim = bwperim(mask);
    perim_thick = imdilate(perim, strel('disk', 1));
    
    cidx = mod(i-1, n_palette) + 1;
    color = base_colors(cidx, :);
    
    for c = 1:3
        ch = white_bg(:,:,c);
        ch(perim_thick) = color(c);
        white_bg(:,:,c) = ch;
    end
end

imshow(white_bg, 'Border', 'tight');
hold on;

% Scale bar on white bg uses black
line([sb_x, sb_x + scale_bar_px], [sb_y, sb_y], 'Color', 'k', 'LineWidth', 4);
text(sb_x + scale_bar_px/2, sb_y - h*0.02, sprintf('%d µm', scale_bar_um), ...
    'Color', 'black', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

text(w * 0.95, h * 0.94, sprintf('n = %d ROIs', nROIs), ...
    'Color', 'black', 'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

title(sprintf('%s: %d ROIs - Outlines (v4.5)', basename, nROIs), ...
    'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

exportgraphics_safe(fig2, fullfile(figDir, [basename '_roi_outlines_v45.png']), 300);
if config.generate_pdfs
    try
        set(fig2, 'PaperPositionMode', 'auto');
        print(fig2, fullfile(figDir, [basename '_roi_outlines_v45.pdf']), '-dpdf', '-bestfit');
    catch
    end
end
close(fig2);

fprintf('    ✓ ROI map + outline map saved\n');

end

%% Helper: Save figure at specified DPI (fallback-safe)
function exportgraphics_safe(fig, filepath, dpi)
try
    exportgraphics(fig, filepath, 'Resolution', dpi);
catch
    saveas(fig, filepath);
end
end

%% =====================================================================
%  v4.4 NEW FIGURE: ROI GALLERY
%  =====================================================================
function generate_roi_gallery_v45(roi_data, figDir, basename, config)
% Generate a gallery of zoomed individual ROI thumbnails
% Shows each ROI with its trace snippet, SNR, and size

nROIs = size(roi_data.roi_centers, 1);
if nROIs == 0
    return;
end

cols = config.fig.gallery_cols;
rows = config.fig.gallery_rows;
rois_per_page = rows * cols;
n_pages = ceil(nROIs / rois_per_page);

% Background image - enhanced contrast
meanImg_masked = roi_data.meanImg .* roi_data.brain_mask;
bg_min = prctile(meanImg_masked(roi_data.brain_mask), 1);
bg_max = prctile(meanImg_masked(roi_data.brain_mask), 99);
meanImg_norm = (meanImg_masked - bg_min) / (bg_max - bg_min + eps);
meanImg_norm = max(0, min(1, meanImg_norm));
meanImg_norm = meanImg_norm .^ 0.45;  % Gamma brighten for thumbnails
[h, w] = size(meanImg_norm);

% Thumbnail size - larger window for context (at least 20px radius)
thumb_half = max(20, round(config.roi.gSig * 6));

for page = 1:n_pages
    fig = figure('Position', [50 50 1600 1000], 'Visible', 'off');
    
    start_roi = (page-1) * rois_per_page + 1;
    end_roi = min(page * rois_per_page, nROIs);
    
    for idx = 1:(end_roi - start_roi + 1)
        roi_num = start_roi + idx - 1;
        
        subplot(rows, cols, idx);
        
        % Get ROI center and crop region
        cx = round(roi_data.roi_centers(roi_num, 1));
        cy = round(roi_data.roi_centers(roi_num, 2));
        
        y1 = max(1, cy - thumb_half);
        y2 = min(h, cy + thumb_half);
        x1 = max(1, cx - thumb_half);
        x2 = min(w, cx + thumb_half);
        
        % Crop image and mask
        thumb_img = meanImg_norm(y1:y2, x1:x2);
        thumb_mask = roi_data.roi_masks(y1:y2, x1:x2, roi_num);
        
        % Build RGB thumbnail with ROI highlight
        thumb_rgb = repmat(thumb_img, [1 1 3]);
        % Tint ROI region green
        perim = bwperim(thumb_mask);
        perim_thick = imdilate(perim, strel('disk', 1));
        for c = 1:3
            ch = thumb_rgb(:,:,c);
            % Subtle green fill
            if c == 2  % Green channel
                ch(thumb_mask) = ch(thumb_mask) * 0.6 + 0.4;
            else
                ch(thumb_mask) = ch(thumb_mask) * 0.7;
            end
            % Bright boundary
            if c == 2
                ch(perim_thick) = 1.0;
            else
                ch(perim_thick) = 0.2;
            end
            thumb_rgb(:,:,c) = ch;
        end
        
        % Display
        imshow(thumb_rgb);
        axis image; axis off; hold on;
        
        % Calculate physical size
        area_px = sum(roi_data.roi_masks(:,:,roi_num), 'all');
        area_um2 = area_px * config.pixelSize_um^2;
        diameter_um = 2 * sqrt(area_um2 / pi);
        
        % Title with info
        title(sprintf('ROI %d (SNR=%.1f, Ø=%.1fµm)', ...
            roi_num, roi_data.snr(roi_num), diameter_um), 'FontSize', 8);
    end
    
    % Page title
    if n_pages > 1
        sgtitle(sprintf('%s: ROI Gallery (%d/%d)', basename, page, n_pages), ...
            'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
        filename = sprintf('%s_roi_gallery_v45_%d.png', basename, page);
    else
        sgtitle(sprintf('%s: ROI Gallery', basename), ...
            'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
        filename = [basename '_roi_gallery_v45.png'];
    end
    
    saveas(fig, fullfile(figDir, filename));
    if config.generate_pdfs
        try
            pdf_filename = strrep(filename, '.png', '.pdf');
            print(fig, fullfile(figDir, pdf_filename), '-dpdf', '-bestfit');
        catch
        end
    end
    close(fig);
end

fprintf('    ✓ ROI gallery saved (%d pages)\n', n_pages);

end

%% =====================================================================
%  TRACE FIGURES (v4.5)
%  =====================================================================
function generate_trace_figures_v45(roi_data, figDir, basename, config)

nROIs = size(roi_data.dFF, 1);
time_vec = (0:size(roi_data.dFF,2)-1) / config.frameRate;

n_per_page = 16;
n_pages = ceil(nROIs / n_per_page);

for page = 1:n_pages
    fig = figure('Position', [50 50 1400 900], 'Visible', 'off');
    
    start_roi = (page-1) * n_per_page + 1;
    end_roi = min(page * n_per_page, nROIs);
    
    for i = 1:(end_roi - start_roi + 1)
        r = start_roi + i - 1;
        subplot(4, 4, i);
        
        if isfield(roi_data, 'dFF_smooth')
            plot(time_vec, roi_data.dFF_smooth(r,:), 'k', 'LineWidth', 0.5);
        else
            plot(time_vec, roi_data.dFF(r,:), 'k', 'LineWidth', 0.5);
        end
        xlabel('Time (s)'); ylabel('dF/F');
        
        % Calculate diameter
        area_px = sum(roi_data.roi_masks(:,:,r), 'all');
        diameter_um = 2 * sqrt(area_px * config.pixelSize_um^2 / pi);
        
        title(sprintf('ROI %d (SNR:%.1f, Ø:%.1fµm)', r, roi_data.snr(r), diameter_um), 'FontSize', 8);
        xlim([0 time_vec(end)]); grid on;
    end
    
    if n_pages > 1
        sgtitle(sprintf('%s: Traces (%d/%d)', basename, page, n_pages), 'Interpreter', 'none');
        pngFile = fullfile(figDir, sprintf('%s_traces_v45_%d.png', basename, page));
    else
        sgtitle(sprintf('%s: Traces', basename), 'Interpreter', 'none');
        pngFile = fullfile(figDir, [basename '_traces_v45.png']);
    end
    
    saveas(fig, pngFile);
    if config.generate_pdfs
        try
            pdfFile = strrep(pngFile, '.png', '.pdf');
            print(fig, pdfFile, '-dpdf', '-bestfit');
        catch
        end
    end
    close(fig);
end

end

%% =====================================================================
%  CSV EXPORT (v4.4 with physical units)
%  =====================================================================
function export_roi_csv_v45(outputDir, basename, roi_data, config)

if size(roi_data.roi_centers, 1) == 0
    % Write diagnostic CSV
    info = roi_data.detection_info;
    pnr_vals = roi_data.PNR(roi_data.brain_mask);
    pnr_prctiles = prctile(pnr_vals, [25, 50, 75, 90, 95]);
    
    T = table();
    T.Metric = {'n_rois'; 'cn_median'; 'cn_range'; 'cn_max'; ...
                'pnr_median'; 'pnr_max'; ...
                'min_pnr_thresh'; 'min_corr_thresh'; 'min_area'; 'max_area'; ...
                'n_seeds'; 'data_quality'; 'gSig'};
    
    T.Value = {'0'; ...
               sprintf('%.4f', info.cn_median); ...
               sprintf('%.4f', info.cn_range); ...
               sprintf('%.4f', info.cn_max); ...
               sprintf('%.2f', median(pnr_vals)); ...
               sprintf('%.2f', info.pnr_max); ...
               sprintf('%.1f', config.roi.min_pnr); ...
               sprintf('%.2f', config.roi.min_corr); ...
               sprintf('%d', config.roi.min_area_px); ...
               sprintf('%d', config.roi.max_area_px); ...
               sprintf('%d', info.n_seeds); ...
               info.data_quality; ...
               sprintf('%d', config.roi.gSig)};
    
    writetable(T, fullfile(outputDir, [basename '_diagnostic_v45.csv']));
    return;
end

% Full ROI metrics CSV with physical units
T = table();
nROIs_out = size(roi_data.roi_masks, 3);

T.ROI_ID = (1:nROIs_out)';
T.Center_X = roi_data.roi_centers(:, 1);
T.Center_Y = roi_data.roi_centers(:, 2);

% Calculate areas and shape metrics with physical units
areas = zeros(nROIs_out, 1);
eccentricity = zeros(nROIs_out, 1);
solidity = zeros(nROIs_out, 1);

for i = 1:nROIs_out
    mask_i = roi_data.roi_masks(:,:,i);
    areas(i) = sum(mask_i, 'all');
    
    props = regionprops(mask_i, 'Eccentricity', 'Solidity');
    if ~isempty(props)
        eccentricity(i) = props(1).Eccentricity;
        solidity(i) = props(1).Solidity;
    end
end

T.Area_px = areas;
T.Area_um2 = areas * config.pixelSize_um^2;  % v4.4: Physical area
T.Diameter_um = 2 * sqrt(T.Area_um2 / pi);   % v4.4: Equivalent diameter
T.Eccentricity = eccentricity;
T.Solidity = solidity;

if isfield(roi_data.detection_info, 'roi_levels')
    T.Detection_Level = roi_data.detection_info.roi_levels;
    level_names = {'strict', 'medium', 'relaxed', 'pnr_only', 'mean_fallback', 'std_fallback'};
    T.Detection_Name = cell(size(T.Detection_Level));
    for i = 1:length(T.Detection_Level)
        lev = T.Detection_Level(i);
        if lev >= 1 && lev <= 6
            T.Detection_Name{i} = level_names{lev};
        else
            T.Detection_Name{i} = 'unknown';
        end
    end
end

T.SNR = roi_data.snr;
T.Event_Rate_Hz = roi_data.event_rate;
T.Event_Amplitude = roi_data.event_amplitude;
T.Event_Count = roi_data.event_count;

T.Mean_dFF = mean(roi_data.dFF, 2);
T.Max_dFF = max(roi_data.dFF, [], 2);
T.Std_dFF = std(roi_data.dFF, 0, 2);

T.P25_dFF = prctile(roi_data.dFF, 25, 2);
T.P50_dFF = prctile(roi_data.dFF, 50, 2);
T.P75_dFF = prctile(roi_data.dFF, 75, 2);
T.P90_dFF = prctile(roi_data.dFF, 90, 2);
T.P95_dFF = prctile(roi_data.dFF, 95, 2);

writetable(T, fullfile(outputDir, [basename '_roi_metrics_v45.csv']));

end

%% =====================================================================
%  CONDITION SUMMARY
%  =====================================================================
function write_condition_summary(results, outputDir, condName, config, versions)

params = struct();
params.condition = condName;
params.pipeline_version = versions.pipeline_version;
params.n_files = height(results);
params.total_rois = sum(results.nROIs);
params.mean_snr = mean(results.meanSNR, 'omitnan');
params.roi_params = config.roi;
params.roi_params.gSig_um = config.roi.gSig * config.pixelSize_um;
params.roi_params.min_area_um2 = config.roi.min_area_px * config.pixelSize_um^2;
params.roi_params.max_area_um2 = config.roi.max_area_px * config.pixelSize_um^2;

jsonFile = fullfile(outputDir, sprintf('roi_params_%s_v45.json', lower(strrep(condName, ' ', '_'))));
fid = fopen(jsonFile, 'w');
if fid > 0
    fwrite(fid, jsonencode(params));
    fclose(fid);
end

end

%% =====================================================================
%  FINAL SUMMARY
%  =====================================================================
function generate_final_summary(results, outputDir, config, versions)

writetable(results, fullfile(outputDir, 'pipeline_summary_v45.csv'));

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                      FINAL SUMMARY (v4.5)                        ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('Detection params: gSig=%d (%.2f µm), min_area=%d, max_area=%d\n', ...
    config.roi.gSig, config.roi.gSig * config.pixelSize_um, config.roi.min_area_px, config.roi.max_area_px);
fprintf('Files: %d | ROIs: %d | Mean SNR: %.2f | Mean Rate: %.3f Hz\n', ...
    height(results), sum(results.nROIs), mean(results.meanSNR,'omitnan'), mean(results.meanEventRate,'omitnan'));

if ismember('condition', results.Properties.VariableNames)
    conditions = unique(string(results.condition));
    for i = 1:length(conditions)
        idx = string(results.condition) == conditions(i);
        fprintf('  %s: %d files, %d ROIs\n', conditions(i), sum(idx), sum(results.nROIs(idx)));
    end
end

fprintf('Results: %s\n', fullfile(outputDir, 'pipeline_summary_v45.csv'));

end

%% =====================================================================
%  UTILITY FUNCTIONS
%  =====================================================================

function config = merge_with_defaults(user_config)
% Merge a user-provided config struct with defaults. Unknown fields
% trigger a warning. Missing fields are filled from defaults.
config = pipeline_config();
config = merge_struct_recursive(config, user_config);
end

function base = merge_struct_recursive(base, override)
% Recursively merge override fields into base struct.
fields = fieldnames(override);
base_fields = fieldnames(base);
for i = 1:numel(fields)
    if ismember(fields{i}, base_fields) && isstruct(base.(fields{i})) && isstruct(override.(fields{i})) ...
            && ~isempty(fieldnames(base.(fields{i})))
        % Special handling: conditions is a struct array, not nested config
        if strcmp(fields{i}, 'conditions')
            base.(fields{i}) = override.(fields{i});
        else
            base.(fields{i}) = merge_struct_recursive(base.(fields{i}), override.(fields{i}));
        end
    elseif ismember(fields{i}, base_fields)
        base.(fields{i}) = override.(fields{i});
    else
        warning('pipeline:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See pipeline_config() for valid fields.', ...
            fields{i});
    end
end
end

function conditions = prompt_for_conditions()
% Prompt the user to define conditions via folder pickers.
if ~usejava('desktop')
    error(['No conditions defined and no GUI available.\n' ...
           'Usage:\n' ...
           '  config = pipeline_config();\n' ...
           '  config.conditions(1).name = ''Control'';\n' ...
           '  config.conditions(1).tifDir = ''/path/to/tif'';\n' ...
           '  config.conditions(1).aviDir = ''/path/to/avi'';\n' ...
           '  COMPLETE_PIPELINE_RUN_ME_v4_5(config);\n\n' ...
           'See also: pipeline_config']);
end

conditions = struct('name', {}, 'tifDir', {}, 'aviDir', {});
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
    
    tifDir = uigetdir('', sprintf('Select TIF folder for %s', name{1}));
    if tifDir == 0, break; end
    
    aviDir = uigetdir('', sprintf('Select AVI/shifts folder for %s', name{1}));
    if aviDir == 0, break; end
    
    conditions(condIdx).name = name{1};
    conditions(condIdx).tifDir = tifDir;
    conditions(condIdx).aviDir = aviDir;
end

if isempty(conditions)
    error('No conditions defined. Exiting.');
end
end

function logf(fid, varargin)
fprintf(varargin{:});
if fid > 0, fprintf(fid, varargin{:}); end
end

function validate_config(config)
% v4.4: Validate configuration parameters to catch common errors early

% Check ROI size parameters
if config.roi.gSig < 1 || config.roi.gSig > 20
    warning('gSig=%d seems unusual. Expected 1-20 for typical neuron imaging.', config.roi.gSig);
end

if config.roi.min_area_px >= config.roi.max_area_px
    error('min_area_px (%d) must be less than max_area_px (%d)', ...
        config.roi.min_area_px, config.roi.max_area_px);
end

% Check physical sizes
min_diameter_um = 2 * sqrt(config.roi.min_area_px * config.pixelSize_um^2 / pi);
max_diameter_um = 2 * sqrt(config.roi.max_area_px * config.pixelSize_um^2 / pi);

if min_diameter_um > 5
    warning('Minimum ROI diameter (%.1f µm) seems large for neurons.', min_diameter_um);
end

if max_diameter_um > 10
    warning('Maximum ROI diameter (%.1f µm) seems large for Drosophila neurons.', max_diameter_um);
end

% Check SNR threshold
if config.roi.min_snr < 1
    warning('min_snr=%.1f is very low. May detect noise as neurons.', config.roi.min_snr);
end

% Check frame rate
if config.frameRate < 1 || config.frameRate > 1000
    warning('frameRate=%.1f Hz seems unusual for calcium imaging.', config.frameRate);
end

fprintf('✓ Configuration validated\n');
end

function results = append_result(results, result)
if isempty(results)
    results = struct2table(result);
else
    results = [results; struct2table(result)];
end
end

function versions = detect_versions()
try
    v = ver('MATLAB');
    if ~isempty(v) && ~isempty(v.Release)
        rel = regexprep(v.Release, '[()]', '');
    else
        rel = version;
    end
catch
    rel = 'unknown';
end

bfVer = 'unknown';
try
    if exist('loci.formats.FormatTools', 'class')
        bfVer = char(loci.formats.FormatTools.VERSION);
    end
catch
end

versions = struct('MATLAB_release', rel, 'BioFormats_version', bfVer, ...
    'pipeline_version', '4.5', 'pipeline_date', '2026-02');
end

function dirs = create_output_dirs(baseDir)
dirs = struct();
dirs.motion = fullfile(baseDir, 'motion_corrected');
dirs.rois = fullfile(baseDir, 'roi_analysis');
dirs.figures = fullfile(baseDir, 'figures');

fields = fieldnames(dirs);
for i = 1:length(fields)
    if ~exist(dirs.(fields{i}), 'dir')
        mkdir(dirs.(fields{i}));
    end
end
end

function files = find_input_files(folder, file_pattern)
if nargin < 2, file_pattern = 'S*.tif'; end
[~, ~, pat_ext] = fileparts(file_pattern);
if isempty(pat_ext)
    % No extension in pattern, search both .tif and .tiff
    d = [dir(fullfile(folder, [file_pattern '.tif'])); dir(fullfile(folder, [file_pattern '.tiff']))];
else
    d = dir(fullfile(folder, file_pattern));
    % Also try alternate extension
    if strcmpi(pat_ext, '.tif')
        alt_pattern = regexprep(file_pattern, '\.tif$', '.tiff', 'ignorecase');
        d = [d; dir(fullfile(folder, alt_pattern))];
    end
end
[~, idx] = unique(lower({d.name}));
d = d(idx);

keep = ~contains({d.name}, '_rigid', 'IgnoreCase', true);
files = d(keep);

if ~isempty(files)
    nums = zeros(length(files), 1);
    for i = 1:length(files)
        tok = regexp(files(i).name, 'S(\d+)', 'tokens');
        if ~isempty(tok), nums(i) = str2double(tok{1}{1}); end
    end
    [~, order] = sort(nums);
    files = files(order);
end
end

function setup_bioformats()
if isempty(which('bfGetReader'))
    % Bio-Formats not on path. Try common install locations.
    % Users should add Bio-Formats to their MATLAB path before running.
    % Download from: https://www.openmicroscopy.org/bio-formats/
    common_paths = {};
    
    % Check user's MATLAB path for bfmatlab directories
    userpath_dir = userpath;
    if ~isempty(userpath_dir)
        candidate = fullfile(userpath_dir, 'bfmatlab');
        if exist(candidate, 'dir')
            common_paths{end+1} = candidate;
        end
    end
    
    for i = 1:length(common_paths)
        if exist(common_paths{i}, 'dir')
            addpath(genpath(common_paths{i}));
            break;
        end
    end
end

if isempty(which('bfGetReader'))
    error(['Bio-Formats not found on MATLAB path.\n' ...
           'Install from: https://www.openmicroscopy.org/bio-formats/\n' ...
           'Then add it: addpath(genpath(''/path/to/bfmatlab''))']);
end

fprintf('✓ Bio-Formats ready\n');
end
