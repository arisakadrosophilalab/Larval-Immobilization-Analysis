function compare_motion_groups(controlDir, experDir, outDir, varargin)
% COMPARE_MOTION_GROUPS  Publication-grade group comparison for GFP movement data.
%
%   Compares two groups of motion-tracked recordings using three core metrics:
%     Mean Speed   — total activity level (distance / time)
%     Median Step  — typical frame-to-frame displacement (robust to outliers)
%     Jitter Ratio — movement pattern quality (1 = straight, high = wandering)
%
%   USAGE:
%     compare_motion_groups()                           % Opens folder picker
%     compare_motion_groups('/ctrl', '/exp', '/out')    % Specify paths
%
%     cfg = compare_motion_config();                    % Get default config
%     cfg.time_window_min = [0 30];                     % First 30 min only
%     cfg.control_label = 'Saline';
%     compare_motion_groups('/ctrl', '/exp', '/out', cfg)
%
%   INPUTS:
%     controlDir - Path to control group folder (contains motion_qc.csv, etc.)
%     experDir   - Path to experimental group folder
%     outDir     - Path for output files (created if needed)
%     cfg        - (optional) Config struct from compare_motion_config()
%
%   OUTPUTS (written to outDir):
%     group_compare_summary.csv   Statistical comparison table
%     per_file_metrics.csv        Per-file metric values
%     sample_rankings_*.csv       Ranked sample table
%     included_files_*.csv        Files that passed QC
%     *.png, *.pdf                Publication-quality figures
%
%   EFFECT SIZE INTERPRETATION (Hedges' g, bias-corrected Cohen's d):
%     |g| < 0.2  = negligible
%     |g| 0.2-0.5 = small
%     |g| 0.5-0.8 = medium
%     |g| >= 0.8  = large
%
%   CONFIDENCE INTERVAL INTERPRETATION:
%     If 95% CI excludes zero: effect is statistically significant at p < 0.05
%     CI width indicates precision: narrower = more precise estimate
%
%   DEPENDENCIES:
%     - MATLAB R2020b or later
%     - Statistics and Machine Learning Toolbox (for ranksum, lillietest, iqr)
%
%   See also: compare_motion_config, analyze_motion_and_QC

%% ==================== PARSE CONFIG ====================
if nargin >= 4 && isstruct(varargin{1})
    cfg = merge_with_defaults_cmp(varargin{1});
elseif nargin >= 4
    error('Fourth argument must be a config struct from compare_motion_config().');
else
    cfg = compare_motion_config();
end

controlLabel  = char(cfg.control_label);
experLabel    = char(cfg.experimental_label);
exportDPI     = cfg.export_dpi;
exportVector  = cfg.export_vector;
robustYLim    = cfg.robust_ylim;
minSampleSize = cfg.min_sample_size;
nB            = cfg.n_bootstrap;
rngSeed       = cfg.rng_seed;
timeWindow    = cfg.time_window_min;

rng(rngSeed, 'twister');

%% ==================== RESOLVE PATHS ====================
if nargin < 1 || isempty(controlDir)
    controlDir = pick_folder('Select CONTROL group folder');
end
if nargin < 2 || isempty(experDir)
    experDir = pick_folder('Select EXPERIMENTAL group folder');
end
if nargin < 3 || isempty(outDir)
    outDir = pick_folder('Select OUTPUT folder');
end
if ~isfolder(outDir), mkdir(outDir); end

%% ==================== HEADER ====================
fprintf('\n');
fprintf('========================================================================\n');
fprintf('     GFP Movement Group Comparison (3 Core Metrics, Pub Quality)        \n');
fprintf('========================================================================\n');
fprintf('\n');
fprintf('%s:\n  %s\n', controlLabel, controlDir);
fprintf('%s:\n  %s\n', experLabel, experDir);
fprintf('Output:\n  %s\n', outDir);
if isinf(timeWindow(2))
    fprintf('Time window: full recording\n\n');
else
    fprintf('Time window: %.0f – %.0f min\n\n', timeWindow(1), timeWindow(2));
end

%% ==================== LOAD DATA ====================
C = load_tables(controlDir);
E = load_tables(experDir);

C = select_pass(C);
E = select_pass(E);

nC = height(C.qc);
nE = height(E.qc);
fprintf('QC PASS: %s = %d, %s = %d\n', controlLabel, nC, experLabel, nE);

%% ==================== SAMPLE SIZE WARNINGS ====================
sampleSizeWarnings = check_sample_sizes(nC, nE, controlLabel, experLabel, minSampleSize);

ctrlSafe = sanitize_filename(controlLabel);
expSafe  = sanitize_filename(experLabel);
if nC > 0
    writetable(C.qc(:,{'File'}), fullfile(outDir, sprintf('included_files_%s.csv', ctrlSafe)));
end
if nE > 0
    writetable(E.qc(:,{'File'}), fullfile(outDir, sprintf('included_files_%s.csv', expSafe)));
end

%% ==================== EXTRACT METRICS ====================
MC = metrics_from_tables(C, controlLabel, controlDir, timeWindow);
ME = metrics_from_tables(E, experLabel, experDir, timeWindow);

%% ==================== DURATION CONSISTENCY CHECK ====================
durationWarnings = check_duration_consistency(MC, ME, controlLabel, experLabel);

%% ==================== STATISTICS (3 CORE METRICS) ====================
S = make_summary_stats(MC, ME, controlLabel, experLabel, nB);

S.AnalysisDate         = repmat(string(datestr(now, 'yyyy-mm-dd HH:MM:SS')), height(S), 1);
S.RNGSeed              = repmat(rngSeed, height(S), 1);
S.BootstrapIterations  = repmat(nB, height(S), 1);
S.TimeWindow_start_min = repmat(timeWindow(1), height(S), 1);
S.TimeWindow_end_min   = repmat(timeWindow(2), height(S), 1);

writetable(S, fullfile(outDir, 'group_compare_summary.csv'));
fprintf('\nWrote: group_compare_summary.csv\n');

allMetrics = [MC; ME];
writetable(allMetrics, fullfile(outDir, 'per_file_metrics.csv'));
fprintf('Wrote: per_file_metrics.csv\n\n');

%% ==================== COLORS ====================
c1 = [0.12 0.47 0.71];   % Blue  (control)
c2 = [1.00 0.50 0.05];   % Orange (experimental)
fnSuffix = sprintf('_%s_vs_%s', expSafe, ctrlSafe);

%% ==================== GENERATE FIGURES ====================
fprintf('Generating figures (3 core metrics)...\n');
mu = char(181);

% --- 3 Core Metric Box-Swarm Plots ---
plot_box_swarm(MC, ME, S, 'Mean_speed_um_per_s', ['Mean speed (' mu 'm/s)'], ...
    fullfile(outDir, ['mean_speed_um_per_s' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector, 'linear', robustYLim);

plot_box_swarm(MC, ME, S, 'Median_step_um', ['Median step (' mu 'm)'], ...
    fullfile(outDir, ['median_step_um' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector, 'linear', robustYLim);

plot_box_swarm(MC, ME, S, 'Jitter_ratio', 'Jitter ratio', ...
    fullfile(outDir, ['jitter_ratio' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector, 'log', false);

% --- Timecourse Overlay ---
plot_timecourse(controlDir, experDir, fullfile(outDir, ['timecourse_overlay' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, nC, nE, exportDPI, exportVector);

% --- Polar Drift Direction ---
plot_polar_drift(MC, ME, fullfile(outDir, ['drift_angle_polar' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector);

% --- Forest Plot ---
plot_forest(S, fullfile(outDir, ['effects_forest' fnSuffix]), ...
    controlLabel, experLabel, exportDPI, exportVector);

% --- Stats Table ---
plot_stats_table(MC, ME, S, fullfile(outDir, ['stats_table' fnSuffix]), ...
    controlLabel, experLabel, exportDPI, exportVector);

% --- Conclusions Panel ---
plot_conclusions(S, fullfile(outDir, ['conclusions_panel' fnSuffix]), ...
    controlLabel, experLabel, exportDPI, exportVector, sampleSizeWarnings, durationWarnings);

% --- Ranked Sample Plots ---
fprintf('Generating ranked sample plots...\n');
plot_ranked_samples(MC, ME, 'Mean_speed_um_per_s', ['Mean speed (' mu 'm/s)'], ...
    fullfile(outDir, ['ranked_mean_speed' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector);

plot_ranked_samples(MC, ME, 'Median_step_um', ['Median step (' mu 'm)'], ...
    fullfile(outDir, ['ranked_median_step' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector);

plot_ranked_samples(MC, ME, 'Jitter_ratio', 'Jitter ratio', ...
    fullfile(outDir, ['ranked_jitter_ratio' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector);

create_ranked_csv(MC, ME, fullfile(outDir, ['sample_rankings' fnSuffix '.csv']), ...
    controlLabel, experLabel);

% --- Drift Rate vs Mean Speed Scatter ---
plot_rate_comparison(MC, ME, fullfile(outDir, ['drift_vs_speed_comparison' fnSuffix]), ...
    c1, c2, controlLabel, experLabel, exportDPI, exportVector);

%% ==================== DONE ====================
fprintf('\n');
fprintf('========================================================================\n');
fprintf('                        ANALYSIS COMPLETE                               \n');
fprintf('========================================================================\n');
fprintf('Output: %d DPI PNG', exportDPI);
if exportVector, fprintf(' + PDF'); end
fprintf('\n');

if ~isempty(sampleSizeWarnings) || ~isempty(durationWarnings)
    fprintf('\n--- WARNINGS ---\n');
    for i = 1:length(sampleSizeWarnings), fprintf('  %s\n', sampleSizeWarnings{i}); end
    for i = 1:length(durationWarnings),   fprintf('  %s\n', durationWarnings{i}); end
end

fprintf('\nOutput: %s\n', outDir);
open_folder(outDir);
end

%% ========================================================================
%% ========================== CONFIG HELPERS ==============================
%% ========================================================================

function folder = pick_folder(prompt)
if usejava('desktop')
    folder = uigetdir('', prompt);
    if folder == 0
        error('No folder selected. Exiting.');
    end
else
    error(['No folder specified and no GUI available.\n' ...
           'Usage: compare_motion_groups(controlDir, experDir, outDir)\n' ...
           'See also: compare_motion_config']);
end
end

function cfg = merge_with_defaults_cmp(user_cfg)
cfg = compare_motion_config();
default_fields = fieldnames(cfg);
user_fields = fieldnames(user_cfg);
for i = 1:numel(user_fields)
    if ismember(user_fields{i}, default_fields)
        cfg.(user_fields{i}) = user_cfg.(user_fields{i});
    else
        warning('compare:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See compare_motion_config().', ...
            user_fields{i});
    end
end
end

%% ========================================================================
%% ========================== SAMPLE SIZE CHECKS ==========================
%% ========================================================================

function warnings = check_sample_sizes(nC, nE, ctrlLabel, expLabel, minN)
warnings = {};

if nC < 3
    warnings{end+1} = sprintf('CRITICAL: %s n=%d (stats unreliable)', ctrlLabel, nC);
    fprintf('*** CRITICAL: %s n=%d - stats unreliable ***\n', ctrlLabel, nC);
end
if nE < 3
    warnings{end+1} = sprintf('CRITICAL: %s n=%d (stats unreliable)', expLabel, nE);
    fprintf('*** CRITICAL: %s n=%d - stats unreliable ***\n', expLabel, nE);
end

if nC >= 3 && nC < minN
    warnings{end+1} = sprintf('CAUTION: %s n=%d (recommend n>=%d)', ctrlLabel, nC, minN);
    fprintf('CAUTION: %s n=%d (recommend >=%d)\n', ctrlLabel, nC, minN);
end
if nE >= 3 && nE < minN
    warnings{end+1} = sprintf('CAUTION: %s n=%d (recommend n>=%d)', expLabel, nE, minN);
    fprintf('CAUTION: %s n=%d (recommend >=%d)\n', expLabel, nE, minN);
end

if nC > 0 && nE > 0
    ratio = max(nC, nE) / min(nC, nE);
    if ratio > 3
        warnings{end+1} = sprintf('IMBALANCE: %d vs %d (ratio=%.1f)', nC, nE, ratio);
    end
end

if isempty(warnings)
    fprintf('Sample sizes OK (n=%d, n=%d)\n', nC, nE);
end
end

%% ========================================================================
%% ====================== DURATION CONSISTENCY CHECK ======================
%% ========================================================================

function warnings = check_duration_consistency(MC, ME, ctrlLabel, expLabel)
warnings = {};

if ~ismember('Duration_s', MC.Properties.VariableNames)
    return;
end

durC = MC.Duration_s(isfinite(MC.Duration_s));
durE = ME.Duration_s(isfinite(ME.Duration_s));

if isempty(durC) && isempty(durE), return; end

if numel(durC) >= 2
    cvC = std(durC) / mean(durC) * 100;
    if cvC > 20
        warnings{end+1} = sprintf('%s: duration CV=%.0f%%', ctrlLabel, cvC);
    else
        fprintf('%s durations: %.0f +/- %.0f s\n', ctrlLabel, mean(durC), std(durC));
    end
end

if numel(durE) >= 2
    cvE = std(durE) / mean(durE) * 100;
    if cvE > 20
        warnings{end+1} = sprintf('%s: duration CV=%.0f%%', expLabel, cvE);
    else
        fprintf('%s durations: %.0f +/- %.0f s\n', expLabel, mean(durE), std(durE));
    end
end

if ~isempty(durC) && ~isempty(durE)
    pctDiff = abs(mean(durC) - mean(durE)) / ((mean(durC) + mean(durE))/2) * 100;
    if pctDiff > 25
        warnings{end+1} = sprintf('Duration mismatch: %.0fs vs %.0fs', mean(durC), mean(durE));
    end
end
end

%% ========================================================================
%% ========================== DATA LOADING ================================
%% ========================================================================

function T = load_tables(folder)
T = struct('qc', table(), 'enr', table());

qcPath = fullfile(folder, 'motion_qc.csv');
if exist(qcPath, 'file')
    T.qc = readtable(qcPath, 'VariableNamingRule', 'preserve');
    T.qc = normalize_file_col(T.qc);
end

enrPath = fullfile(folder, 'motion_enriched.csv');
if exist(enrPath, 'file')
    T.enr = readtable(enrPath, 'VariableNamingRule', 'preserve');
    T.enr = normalize_file_col(T.enr);
end
end

function T = normalize_file_col(T)
if isempty(T), return; end

if ~ismember('File', T.Properties.VariableNames)
    vars = T.Properties.VariableNames;
    hit = find(strcmpi(vars, 'file') | contains(lower(vars), 'file'), 1);
    if ~isempty(hit)
        T.File = T.(vars{hit});
    else
        T.File = repmat("unknown", height(T), 1);
    end
end

T.File = string(T.File);
T.File = regexprep(T.File, '\.csv$', '.avi', 'ignorecase');
end

function T = select_pass(T)
if isempty(T.qc), return; end

if ismember('QC_PASS', T.qc.Properties.VariableNames)
    passIdx = T.qc.QC_PASS == true | T.qc.QC_PASS == 1;
    T.qc = T.qc(passIdx, :);
else
    warning('QC_PASS column not found; including all rows.');
end

if ~isempty(T.enr)
    T.enr = T.enr(ismember(T.enr.File, T.qc.File), :);
end
end

%% ========================================================================
%% ======================= METRIC EXTRACTION ==============================
%% ========================================================================

function M = metrics_from_tables(T, groupName, dataDir, timeWindow)
% Extract motion metrics from QC + enriched tables.
% If timeWindow is not [0 Inf], recompute metrics from per-frame shift CSVs.

if isempty(T.qc)
    M = table();
    return;
end

qc  = T.qc;
enr = T.enr;
n   = height(qc);

useTimeWindow = ~(timeWindow(1) == 0 && isinf(timeWindow(2)));

M       = table();
M.File  = qc.File;
M.Group = repmat(string(groupName), n, 1);

if useTimeWindow
    % Recompute metrics from per-frame shift CSVs within the time window
    fprintf('  %s: recomputing metrics for [%.0f, %.0f] min window...\n', ...
        groupName, timeWindow(1), timeWindow(2));

    M.Distance_traveled_um = nan(n, 1);
    M.Net_displacement_um  = nan(n, 1);
    M.Drift_rate_um_per_s  = nan(n, 1);
    M.Drift_angle_deg      = nan(n, 1);
    M.Median_step_um       = nan(n, 1);
    M.Jitter_ratio         = nan(n, 1);
    M.Duration_s           = nan(n, 1);
    M.Mean_speed_um_per_s  = nan(n, 1);

    for i = 1:n
        fname = char(qc.File(i));
        stem  = regexprep(fname, '\.[^.]+$', '');
        csvPath = fullfile(dataDir, [stem '_shifts.csv']);
        if ~exist(csvPath, 'file'), continue; end

        Tcsv = readtable(csvPath);

        % Determine time column
        if ismember('time_min', Tcsv.Properties.VariableNames)
            t_min = Tcsv.time_min;
        elseif ismember('time_s', Tcsv.Properties.VariableNames)
            t_min = Tcsv.time_s / 60;
        else
            continue;
        end

        % Apply time window
        mask = t_min >= timeWindow(1) & t_min <= timeWindow(2);
        if sum(mask) < 10, continue; end

        % Extract shifts in the window
        if all(ismember({'dx_um','dy_um'}, Tcsv.Properties.VariableNames))
            dx = Tcsv.dx_um(mask);
            dy = Tcsv.dy_um(mask);
        elseif all(ismember({'dx_px','dy_px','px_um'}, Tcsv.Properties.VariableNames))
            px_um = Tcsv.px_um(1);
            dx = Tcsv.dx_px(mask) * px_um;
            dy = Tcsv.dy_px(mask) * px_um;
        else
            continue;
        end

        % Determine fps
        if ismember('fps_used', Tcsv.Properties.VariableNames)
            fps = Tcsv.fps_used(1);
        else
            fps = 1;
        end

        steps = hypot([0; diff(dx)], [0; diff(dy)]);
        nF = numel(dx);
        dur_s = (nF - 1) / max(1, fps);

        % Median filter for net drift
        win = max(3, 2*floor((60*fps)/2)+1);
        try
            dx_s = medfilt1(dx, win, 'omitnan', 'truncate');
            dy_s = medfilt1(dy, win, 'omitnan', 'truncate');
        catch
            dx_s = movmedian(dx, win, 'omitnan');
            dy_s = movmedian(dy, win, 'omitnan');
        end
        net_drift = hypot(dx_s(end)-dx_s(1), dy_s(end)-dy_s(1));

        cum_path = sum(steps, 'omitnan');

        M.Distance_traveled_um(i) = cum_path;
        M.Net_displacement_um(i)  = net_drift;
        M.Drift_rate_um_per_s(i)  = net_drift / max(1e-9, dur_s);
        M.Drift_angle_deg(i)      = atan2d(dy_s(end)-dy_s(1), dx_s(end)-dx_s(1));
        M.Median_step_um(i)       = median(steps, 'omitnan');
        M.Jitter_ratio(i)         = cum_path / max(1e-9, net_drift);
        M.Duration_s(i)           = dur_s;
        M.Mean_speed_um_per_s(i)  = cum_path / max(1e-9, dur_s);
    end
else
    % Use pre-computed values from the summary tables (original behavior)
    M.Distance_traveled_um = get_col(qc, enr, {'CumPath_postSettle_um', 'CumPath_um'});
    M.Net_displacement_um  = get_col(qc, enr, {'NetDrift_postSettle_um', 'NetDrift_um'});
    M.Drift_rate_um_per_s  = get_col(qc, enr, {'DriftRate_postSettle_um_per_s', 'DriftRate_um_per_s'});
    M.Drift_angle_deg      = get_col(qc, enr, {'DriftAngle_deg'});

    M.Median_step_um = get_col(qc, enr, {'MedianDisp_postSettle_um', 'MedianDisp_um'});

    jit = get_col(qc, enr, {'JitterRatio_post_cumPath_over_netDrift', 'JitterRatio'});
    if all(isnan(jit))
        cp = M.Distance_traveled_um;
        nd = M.Net_displacement_um;
        jit = cp ./ max(nd, 1e-12);
    end
    M.Jitter_ratio = jit;

    dur = get_col(qc, enr, {'Duration_postSettle_s', 'Duration_s', 'RecordingDuration_s', ...
                             'TotalTime_s', 'duration_s', 'duration', 'Time_s'});
    if all(isnan(dur))
        ok = isfinite(M.Net_displacement_um) & isfinite(M.Drift_rate_um_per_s) & M.Drift_rate_um_per_s > 1e-12;
        dur = nan(n, 1);
        dur(ok) = M.Net_displacement_um(ok) ./ M.Drift_rate_um_per_s(ok);
    end
    M.Duration_s = dur;

    M.Mean_speed_um_per_s = nan(n, 1);
    ok = isfinite(M.Jitter_ratio) & isfinite(M.Drift_rate_um_per_s);
    if any(ok)
        M.Mean_speed_um_per_s(ok) = M.Jitter_ratio(ok) .* M.Drift_rate_um_per_s(ok);
    end
end

nValid = sum(isfinite(M.Mean_speed_um_per_s));
fprintf('  %s: %d/%d samples have Mean Speed\n', groupName, nValid, n);
end

function v = get_col(qc, enr, names)
n = height(qc);
v = nan(n, 1);

for i = 1:numel(names)
    nm = names{i};
    if ismember(nm, qc.Properties.VariableNames)
        v = safe_to_double(qc.(nm));
        return;
    end
    if ~isempty(enr) && ismember(nm, enr.Properties.VariableNames)
        [tf, idx] = ismember(qc.File, enr.File);
        vals = safe_to_double(enr.(nm));
        v(tf) = vals(idx(tf));
        return;
    end
end
end

function v = safe_to_double(x)
if isnumeric(x)
    v = double(x);
elseif iscell(x)
    try
        v = cellfun(@(c) str2double(string(c)), x);
    catch
        v = nan(size(x));
    end
elseif isstring(x) || ischar(x)
    v = str2double(x);
else
    v = nan(size(x, 1), 1);
end
v = v(:);
end

%% ========================================================================
%% ========================= STATISTICS ==================================
%% ========================================================================

function S = make_summary_stats(MC, ME, ctrlLabel, expLabel, nB)
metrics = {'Mean_speed_um_per_s', 'Median_step_um', 'Jitter_ratio'};
pretty  = {'Mean speed (um/s)',   'Median step (um)', 'Jitter ratio'};

nM = numel(metrics);

MetricVar               = metrics(:);
Metric                  = pretty(:);
N_Control               = nan(nM, 1);
Mean_Control            = nan(nM, 1);
Median_Control          = nan(nM, 1);
SD_Control              = nan(nM, 1);
SEM_Control             = nan(nM, 1);
IQR_Control             = nan(nM, 1);
N_Experimental          = nan(nM, 1);
Mean_Experimental       = nan(nM, 1);
Median_Experimental     = nan(nM, 1);
SD_Experimental         = nan(nM, 1);
SEM_Experimental        = nan(nM, 1);
IQR_Experimental        = nan(nM, 1);
MeanDiff_ExpMinusCtrl   = nan(nM, 1);
MedianDiff_ExpMinusCtrl = nan(nM, 1);
BootCI_lo               = nan(nM, 1);
BootCI_hi               = nan(nM, 1);
P_Welch                 = nan(nM, 1);
P_MannWhitney           = nan(nM, 1);
P_Permutation           = nan(nM, 1);
Hedges_g                = nan(nM, 1);
Normal_Control          = nan(nM, 1);
Normal_Experimental     = nan(nM, 1);

for i = 1:nM
    if ~ismember(metrics{i}, MC.Properties.VariableNames) || ...
       ~ismember(metrics{i}, ME.Properties.VariableNames)
        continue;
    end

    xc = MC.(metrics{i}); xc = xc(isfinite(xc));
    xe = ME.(metrics{i}); xe = xe(isfinite(xe));
    nc = numel(xc); ne = numel(xe);

    N_Control(i) = nc;
    N_Experimental(i) = ne;

    if nc > 0
        Mean_Control(i)   = mean(xc);
        Median_Control(i) = median(xc);
        SD_Control(i)     = std(xc);
        SEM_Control(i)    = std(xc) / sqrt(nc);
        IQR_Control(i)    = iqr(xc);
    end
    if ne > 0
        Mean_Experimental(i)   = mean(xe);
        Median_Experimental(i) = median(xe);
        SD_Experimental(i)     = std(xe);
        SEM_Experimental(i)    = std(xe) / sqrt(ne);
        IQR_Experimental(i)    = iqr(xe);
    end

    % Normality (Lilliefors)
    if nc >= 4
        try
            warning('off', 'all');
            [~, pN] = lillietest(xc);
            Normal_Control(i) = pN > 0.05;
            warning('on', 'all');
        catch
            Normal_Control(i) = NaN;
        end
    end
    if ne >= 4
        try
            warning('off', 'all');
            [~, pN] = lillietest(xe);
            Normal_Experimental(i) = pN > 0.05;
            warning('on', 'all');
        catch
            Normal_Experimental(i) = NaN;
        end
    end

    if nc < 2 || ne < 2, continue; end

    MeanDiff_ExpMinusCtrl(i)   = mean(xe) - mean(xc);
    MedianDiff_ExpMinusCtrl(i) = median(xe) - median(xc);

    % Welch t-test
    P_Welch(i) = welch_ttest(xc, xe);

    % Mann-Whitney U
    try
        P_MannWhitney(i) = ranksum(xc, xe);
    catch
        P_MannWhitney(i) = NaN;
    end

    % Bootstrap CI for mean difference
    diffs = nan(nB, 1);
    for b = 1:nB
        diffs(b) = mean(xe(randi(ne, ne, 1))) - mean(xc(randi(nc, nc, 1)));
    end
    BootCI_lo(i) = quantile(diffs, 0.025);
    BootCI_hi(i) = quantile(diffs, 0.975);

    % Permutation test
    pool     = [xc; xe];
    diff_obs = abs(mean(xe) - mean(xc));
    count    = 0;
    for b = 1:nB
        idx = randperm(numel(pool));
        d   = abs(mean(pool(idx(1:ne))) - mean(pool(idx(ne+1:end))));
        if d >= diff_obs, count = count + 1; end
    end
    P_Permutation(i) = (count + 1) / (nB + 1);

    % Hedges' g (bias-corrected Cohen's d)
    sp = sqrt(((nc-1)*var(xc) + (ne-1)*var(xe)) / (nc + ne - 2));
    if sp > 0
        d = (mean(xe) - mean(xc)) / sp;
        J = 1 - 3 / (4*(nc + ne) - 9);
        Hedges_g(i) = J * d;
    end
end

P_FDR = fdr_bh(P_Permutation);

S = table(MetricVar, Metric, ...
    N_Control, Mean_Control, Median_Control, SD_Control, SEM_Control, IQR_Control, Normal_Control, ...
    N_Experimental, Mean_Experimental, Median_Experimental, SD_Experimental, SEM_Experimental, IQR_Experimental, Normal_Experimental, ...
    MeanDiff_ExpMinusCtrl, MedianDiff_ExpMinusCtrl, BootCI_lo, BootCI_hi, ...
    P_Welch, P_MannWhitney, P_Permutation, P_FDR, Hedges_g);
S.ControlLabel = repmat(string(ctrlLabel), nM, 1);
S.ExperLabel   = repmat(string(expLabel),  nM, 1);
end

function p = welch_ttest(x, y)
mx = mean(x); my = mean(y);
vx = var(x);  vy = var(y);
nx = numel(x); ny = numel(y);
se2 = vx/nx + vy/ny;
if se2 <= 0, p = 1; return; end
t  = (my - mx) / sqrt(se2);
df = (se2^2) / ((vx^2)/((nx^2)*(nx-1) + eps) + (vy^2)/((ny^2)*(ny-1) + eps) + eps);
p  = 2 * tcdf(-abs(t), max(df, 1));
end

function q = fdr_bh(p)
p = p(:);
m = numel(p);
nanMask = isnan(p);
p(nanMask) = 1;
[ps, idx] = sort(p);
q = nan(m, 1);
for k = m:-1:1
    if k == m
        q(idx(k)) = ps(k);
    else
        q(idx(k)) = min(q(idx(k+1)), ps(k) * m / k);
    end
end
q = min(q, 1);
q(nanMask) = NaN;
end
%% ========================================================================
%% ============================== PLOTTING ================================
%% ========================================================================

function plot_box_swarm(MC, ME, S, field, ylab, pathBase, c1, c2, ctrlLabel, expLabel, dpi, vec, scaleMode, robustY)
if ~ismember(field, MC.Properties.VariableNames) || ~ismember(field, ME.Properties.VariableNames)
    fprintf('  Skipped: %s (not found)\n', field);
    return;
end

xc = MC.(field); xc = xc(isfinite(xc));
xe = ME.(field); xe = xe(isfinite(xe));
nc = numel(xc); ne = numel(xe);

if nc == 0 || ne == 0
    fprintf('  Skipped: %s (empty group)\n', field);
    return;
end

useLog = strcmpi(scaleMode, 'log');

if useLog
    allData = [xc; xe];
    minPos  = min(allData(allData > 0));
    if isempty(minPos), minPos = 1e-6; end
    epsVal  = minPos * 0.1;
    nZerosFixed = sum(xc <= 0) + sum(xe <= 0);
    xc(xc <= 0) = epsVal;
    xe(xe <= 0) = epsVal;
else
    nZerosFixed = 0;
end

fig = figure('Color', 'w', 'Position', [100 100 700 420], 'Visible', 'off');

allData = [xc(:); xe(:)];
grp     = [ones(nc, 1); 2*ones(ne, 1)];

boxplot(allData, grp, 'Labels', {ctrlLabel, expLabel}, 'Symbol', '', 'Widths', 0.45);
hold on;

h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), [0.85 0.85 0.85], ...
        'FaceAlpha', 0.5, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.2);
end

jw = 0.12;
scatter(1 + jw*(rand(nc,1)-0.5), xc, 55, c1, 'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.6);
scatter(2 + jw*(rand(ne,1)-0.5), xe, 55, c2, 'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.6);

if useLog
    set(gca, 'YScale', 'log');
    yl = ylim;
    ylim([yl(1)*0.8, yl(2)*1.4]);
end

outlierNote = '';
if robustY && ~useLog && ~isempty(allData)
    q1  = prctile(allData, 25);
    q3  = prctile(allData, 75);
    iqr_val    = q3 - q1;
    upperFence = q3 + 1.5 * iqr_val;
    nOut = sum(allData > upperFence);
    if nOut > 0
        ylim([0, upperFence * 1.25]);
        outlierNote = sprintf('%d outlier%s clipped', nOut, ternary(nOut > 1, 's', ''));
    end
end

% Significance annotation
sigStar  = '';
pValNote = '';
metricIdx = find(strcmp(S.MetricVar, field), 1);
if ~isempty(metricIdx)
    pPerm = S.P_Permutation(metricIdx);
    pFDR  = S.P_FDR(metricIdx);
    if isfinite(pFDR)
        if     pFDR < 0.001, sigStar = '***';
        elseif pFDR < 0.01,  sigStar = '**';
        elseif pFDR < 0.05,  sigStar = '*';
        end
        pValNote = sprintf('p = %s, q_{FDR} = %s', p_fmt(pPerm), p_fmt(pFDR));
    end
end

if ~isempty(sigStar)
    yl = ylim;
    if useLog
        starY    = 10^(log10(yl(2)) * 0.90);
        bracketY = 10^(log10(yl(2)) * 0.85);
    else
        starY    = yl(2) * 0.94;
        bracketY = yl(2) * 0.88;
    end
    line([1 2], [bracketY bracketY], 'Color', 'k', 'LineWidth', 1.5);
    line([1 1], [bracketY bracketY*0.97], 'Color', 'k', 'LineWidth', 1.5);
    line([2 2], [bracketY bracketY*0.97], 'Color', 'k', 'LineWidth', 1.5);
    text(1.5, starY, sigStar, 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

ylabel(ylab, 'FontSize', 12, 'FontWeight', 'bold');
title(ylab,  'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'FontName', 'Arial', 'Box', 'off', 'LineWidth', 1.2, ...
    'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabelRotation', 0);

h1 = plot(nan, nan, 'o', 'MarkerFaceColor', c1, 'MarkerEdgeColor', [0.15 0.15 0.15], 'MarkerSize', 8);
h2 = plot(nan, nan, 'o', 'MarkerFaceColor', c2, 'MarkerEdgeColor', [0.15 0.15 0.15], 'MarkerSize', 8);
leg = legend([h1 h2], {sprintf('%s (n=%d)', ctrlLabel, nc), sprintf('%s (n=%d)', expLabel, ne)}, ...
    'Location', 'northeastoutside', 'FontSize', 9);
set(leg, 'Box', 'off');
legPos = get(leg, 'Position');
set(leg, 'Position', [legPos(1)-0.03, legPos(2), legPos(3), legPos(4)]);

noteY = 0.02;
if useLog && nZerosFixed > 0
    text(0.02, noteY, sprintf('Log scale; %d zeros adjusted', nZerosFixed), ...
        'Units', 'normalized', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
    noteY = noteY + 0.04;
end
if ~isempty(outlierNote)
    text(0.02, noteY, outlierNote, 'Units', 'normalized', 'FontSize', 8, 'Color', [0.5 0.3 0.3]);
    noteY = noteY + 0.04;
end
if ~isempty(pValNote)
    text(0.02, noteY, pValNote, 'Units', 'normalized', 'FontSize', 9, 'Color', [0.2 0.2 0.5], 'FontWeight', 'bold');
end

hold off;
save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ========================= TIMECOURSE OVERLAY ===========================
%% ========================================================================

function plot_timecourse(ctrlDir, expDir, pathBase, c1, c2, ctrlLabel, expLabel, nC, nE, dpi, vec)
tcC = load_timecourse(ctrlDir);
tcE = load_timecourse(expDir);

if isempty(tcC) && isempty(tcE)
    fprintf('  Skipped: timecourse (no data)\n');
    return;
end

mu = char(181);
pm = char(177);

fig = figure('Color', 'w', 'Position', [100 100 680 420], 'Visible', 'off');
hold on;

maxN_C = nC; maxN_E = nE;

if ~isempty(tcC)
    x = tcC.minute;
    y = tcC.mean_perMinuteMedian_um;
    e = tcC.sem_perMinuteMedian_um;
    fill([x; flipud(x)], [y-e; flipud(y+e)], c1, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if ismember('n_files_active', tcC.Properties.VariableNames)
        maxN_C = max(tcC.n_files_active);
    end
    plot(x, y, '-', 'Color', c1, 'LineWidth', 2.5, 'DisplayName', sprintf('%s (n=%d)', ctrlLabel, maxN_C));
end

if ~isempty(tcE)
    x = tcE.minute;
    y = tcE.mean_perMinuteMedian_um;
    e = tcE.sem_perMinuteMedian_um;
    fill([x; flipud(x)], [y-e; flipud(y+e)], c2, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    if ismember('n_files_active', tcE.Properties.VariableNames)
        maxN_E = max(tcE.n_files_active);
    end
    plot(x, y, '-', 'Color', c2, 'LineWidth', 2.5, 'DisplayName', sprintf('%s (n=%d)', expLabel, maxN_E));
end

xlabel('Time (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(['Per-minute median step (' mu 'm)'], 'FontSize', 12, 'FontWeight', 'bold');
title(['Motion over time (mean ' pm ' SEM)'], 'FontSize', 13, 'FontWeight', 'bold');
leg = legend('Location', 'northeast', 'FontSize', 10);
set(leg, 'Box', 'off');
set(gca, 'FontSize', 11, 'FontName', 'Arial', 'Box', 'off', 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

function tc = load_timecourse(folder)
tc = [];
tcPath = fullfile(folder, 'group_timecourse.csv');
if ~exist(tcPath, 'file'), return; end

T   = readtable(tcPath, 'VariableNamingRule', 'preserve');
req = {'minute', 'mean_perMinuteMedian_um', 'sem_perMinuteMedian_um'};
if all(ismember(req, T.Properties.VariableNames))
    tc = T;
end
end

%% ========================================================================
%% ========================= POLAR DRIFT ==================================
%% ========================================================================

function plot_polar_drift(MC, ME, pathBase, c1, c2, ctrlLabel, expLabel, dpi, vec)
thC = MC.Drift_angle_deg; thC = thC(isfinite(thC));
thE = ME.Drift_angle_deg; thE = thE(isfinite(thE));
nc = numel(thC); ne = numel(thE);

if nc == 0 && ne == 0
    fprintf('  Skipped: polar drift (no data)\n');
    return;
end

fig = figure('Color', 'w', 'Position', [100 100 650 480], 'Visible', 'off');
ax  = polaraxes;
hold(ax, 'on');

toRad = @(d) d * pi / 180;

if nc > 0
    polarhistogram(ax, toRad(thC), 12, 'Normalization', 'probability', ...
        'FaceColor', c1, 'FaceAlpha', 0.6, 'EdgeColor', c1*0.6, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s (n=%d)', ctrlLabel, nc));
end
if ne > 0
    polarhistogram(ax, toRad(thE), 12, 'Normalization', 'probability', ...
        'FaceColor', c2, 'FaceAlpha', 0.6, 'EdgeColor', c2*0.6, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s (n=%d)', expLabel, ne));
end

leg = legend(ax, 'Location', 'northeastoutside', 'FontSize', 10);
set(leg, 'Box', 'off');
title(ax, 'Net drift direction', 'FontSize', 13, 'FontWeight', 'bold');
ax.FontSize = 10;
ax.FontName = 'Arial';

save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ========================= FOREST PLOT ==================================
%% ========================================================================

function plot_forest(S, pathBase, ctrlLabel, expLabel, dpi, vec)
k = height(S);
if k == 0, return; end

[~, ord] = sort(abs(S.MeanDiff_ExpMinusCtrl), 'descend');
S = S(ord, :);

x    = S.MeanDiff_ExpMinusCtrl;
lo   = S.BootCI_lo;
hi   = S.BootCI_hi;
g    = S.Hedges_g;
pFDR = S.P_FDR;
y    = (1:k)';

fig = figure('Color', 'w', 'Position', [100 100 820 300], 'Visible', 'off');
hold on;

for i = 1:k
    if isfinite(lo(i)) && isfinite(hi(i))
        lineCol = [0.4 0.4 0.4];
        if pFDR(i) < 0.05, lineCol = [0.1 0.6 0.2]; end
        line([lo(i) hi(i)], [y(i) y(i)], 'Color', lineCol, 'LineWidth', 2.5);
    end
end

for i = 1:k
    if pFDR(i) < 0.05
        scatter(x(i), y(i), 130, [0.1 0.6 0.2], 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    else
        scatter(x(i), y(i), 110, [0.35 0.35 0.35], 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    end
end

xline(0, 'k--', 'LineWidth', 1.5);

allX = [x; lo; hi];
allX = allX(isfinite(allX));
if isempty(allX), allX = [-1 1]; end
xr = max(allX) - min(allX);
if xr == 0, xr = 1; end
xlims = [min(allX) - 0.12*xr, max(allX) + 0.35*xr];
xlim(xlims);

set(gca, 'YTick', y, 'YTickLabel', cellstr(S.Metric), 'FontSize', 11, 'FontName', 'Arial');
ylim([0.3, k + 0.7]);

for i = 1:k
    if isfinite(g(i))
        labelTxt = sprintf('g = %.2f,  q = %s', g(i), p_fmt(pFDR(i)));
        if pFDR(i) < 0.05
            labelTxt = [labelTxt '  *'];
            txtCol = [0.1 0.6 0.2];
        else
            txtCol = [0.4 0.4 0.4];
        end
        text(xlims(2) - 0.02*xr, y(i), labelTxt, 'FontSize', 10, 'Color', txtCol, 'HorizontalAlignment', 'right');
    end
end

xlabel(sprintf('Effect (%s  -  %s)', expLabel, ctrlLabel), 'FontSize', 12, 'FontWeight', 'bold');
title('Effect sizes with 95% bootstrap CI', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'Box', 'off', 'LineWidth', 1.2, 'TickDir', 'out');

text(0.02, 0.04, 'Green = q_{FDR} < 0.05', 'Units', 'normalized', 'FontSize', 9, 'Color', [0.1 0.6 0.2]);

hold off;
save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ========================== STATS TABLE =================================
%% ========================================================================

function plot_stats_table(MC, ME, S, pathBase, ctrlLabel, expLabel, dpi, vec)
mu      = char(181);
metrics = {'Mean_speed_um_per_s', 'Median_step_um', 'Jitter_ratio'};
pretty  = {'Mean speed',          'Median step',     'Jitter ratio'};
units   = {[mu '/s'],             mu,                '-'};

fig = figure('Color', 'w', 'Position', [100 100 1050 320], 'Visible', 'off');
axis off;

text(0.5, 0.95, sprintf('Descriptive Statistics:  %s  vs  %s', expLabel, ctrlLabel), ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Units', 'normalized');

cols = [0.02, 0.18, 0.35, 0.52, 0.69, 0.86];
fs = 10;
y  = 0.78;

text(cols(1), y, 'Metric',                             'FontWeight', 'bold', 'FontSize', fs);
text(cols(2), y, sprintf('%s (n)', ctrlLabel),          'FontWeight', 'bold', 'FontSize', fs);
text(cols(3), y, 'Mean / Median / SD',                 'FontWeight', 'bold', 'FontSize', fs);
text(cols(4), y, sprintf('%s (n)', expLabel),           'FontWeight', 'bold', 'FontSize', fs);
text(cols(5), y, 'Mean / Median / SD',                 'FontWeight', 'bold', 'FontSize', fs);
text(cols(6), y, 'p (q_{FDR})',                         'FontWeight', 'bold', 'FontSize', fs);

line([0.02 0.98], [y-0.04 y-0.04], 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

yy = y - 0.16;
for i = 1:numel(metrics)
    if ~ismember(metrics{i}, MC.Properties.VariableNames), continue; end

    xc = MC.(metrics{i}); xc = xc(isfinite(xc));
    xe = ME.(metrics{i}); xe = xe(isfinite(xe));
    if isempty(xc) && isempty(xe), continue; end

    muC = mean(xc); mdC = median(xc); sdC = std(xc);
    muE = mean(xe); mdE = median(xe); sdE = std(xe);

    metricIdx = find(strcmp(S.MetricVar, metrics{i}), 1);
    if ~isempty(metricIdx) && isfinite(S.P_Permutation(metricIdx))
        pPerm = S.P_Permutation(metricIdx);
        pFDR  = S.P_FDR(metricIdx);
        pStr  = sprintf('%s (%s)', p_fmt(pPerm), p_fmt(pFDR));
        if pFDR < 0.05
            pCol = [0.1 0.6 0.2];
            pStr = [pStr ' *'];
        else
            pCol = [0.3 0.3 0.3];
        end
    else
        pStr = 'N/A';
        pCol = [0.5 0.5 0.5];
    end

    metricStr = sprintf('%s (%sm)', pretty{i}, units{i});
    text(cols(1), yy, metricStr, 'FontSize', fs-1);
    text(cols(2), yy, sprintf('%d', numel(xc)), 'FontSize', fs-1);
    text(cols(3), yy, sprintf('%s / %s / %s', smart_fmt(muC), smart_fmt(mdC), smart_fmt(sdC)), 'FontSize', fs-1);
    text(cols(4), yy, sprintf('%d', numel(xe)), 'FontSize', fs-1);
    text(cols(5), yy, sprintf('%s / %s / %s', smart_fmt(muE), smart_fmt(mdE), smart_fmt(sdE)), 'FontSize', fs-1);
    text(cols(6), yy, pStr, 'FontSize', fs-1, 'Color', pCol, 'FontWeight', 'bold');

    yy = yy - 0.16;
end

text(0.02, 0.08, 'p = permutation test (B = 10,000);  q_{FDR} = Benjamini-Hochberg corrected;  * = q < 0.05', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(0.02, 0.02, 'Mean Speed = distance/time (overall activity);  Jitter Ratio = distance/displacement (1 = straight line)', ...
    'FontSize', 9, 'Color', [0.3 0.3 0.6]);

save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ======================== CONCLUSIONS PANEL =============================
%% ========================================================================

function plot_conclusions(S, pathBase, ctrlLabel, expLabel, dpi, vec, sampleWarnings, durationWarnings)
nM = height(S);
if nM == 0, return; end

figW = 1100;
figH = 900;   % Shorter now with only 3 metrics

fig = figure('Color', 'w', 'Position', [50 50 figW figH], 'Visible', 'off');

ax = axes('Units', 'normalized', 'Position', [0.04 0.03 0.92 0.94], 'Visible', 'off');
set(ax, 'XLim', [0 100], 'YLim', [0 120]);
hold on;

leftMargin  = 2;
currentY    = 115;
lineSpacing = 2.6;
metricGap   = 12;

% ==================== TITLE ====================
text(50, currentY, sprintf('Statistical Conclusions:  %s  vs  %s', expLabel, ctrlLabel), ...
    'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
currentY = currentY - 6;

% ==================== WARNINGS ====================
if ~isempty(sampleWarnings) || ~isempty(durationWarnings)
    text(leftMargin, currentY, 'DATA QUALITY NOTES:', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0.4 0]);
    currentY = currentY - lineSpacing;
    allWarnings = [sampleWarnings(:); durationWarnings(:)];
    for w = 1:min(3, numel(allWarnings))
        text(leftMargin + 1, currentY, ['- ' allWarnings{w}], 'FontSize', 9, 'Color', [0.6 0.3 0]);
        currentY = currentY - lineSpacing;
    end
    currentY = currentY - 1;
end

% Track summary
nSig   = 0;
nLarge = 0;
nMed   = 0;
nSmall = 0;
pctReductions = [];

isLowerBetter = @(name) contains(lower(name), ["speed", "step", "jitter", "drift", "distance", "displacement"]);

% ==================== METRICS LOOP ====================
for i = 1:nM
    metric   = S.Metric{i};
    d        = S.MeanDiff_ExpMinusCtrl(i);
    pPerm    = S.P_Permutation(i);
    pFDR     = S.P_FDR(i);
    hg       = S.Hedges_g(i);
    ci_lo    = S.BootCI_lo(i);
    ci_hi    = S.BootCI_hi(i);
    nCtrl    = S.N_Control(i);
    nExp     = S.N_Experimental(i);
    meanCtrl = S.Mean_Control(i);
    meanExp  = S.Mean_Experimental(i);

    pMW = NaN;
    if ismember('P_MannWhitney', S.Properties.VariableNames)
        pMW = S.P_MannWhitney(i);
    end

    % Tracking
    if isfinite(pFDR) && pFDR < 0.05, nSig = nSig + 1; end
    if isfinite(hg)
        if     abs(hg) >= 0.8, nLarge = nLarge + 1;
        elseif abs(hg) >= 0.5, nMed   = nMed + 1;
        elseif abs(hg) >= 0.2, nSmall = nSmall + 1;
        end
    end

    pctChange = NaN;
    if isfinite(meanCtrl) && isfinite(meanExp) && abs(meanCtrl) > 1e-12
        pctChange = ((meanExp - meanCtrl) / abs(meanCtrl)) * 100;
        if isLowerBetter(metric) && pctChange < 0
            pctReductions(end+1) = abs(pctChange); %#ok<AGROW>
        end
    end

    % Color
    if isfinite(pFDR) && pFDR < 0.05
        metricCol = [0.1 0.55 0.2];
        sigMark   = ' *';
    else
        metricCol = [0.15 0.15 0.15];
        sigMark   = '';
    end

    % Header
    text(leftMargin, currentY, sprintf('%s%s   [n = %d vs %d]', metric, sigMark, nCtrl, nExp), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', metricCol);
    currentY = currentY - lineSpacing;

    if isfinite(d) && isfinite(ci_lo) && isfinite(ci_hi) && isfinite(hg)

        % Line 1: Mean diff + CI
        text(leftMargin + 2, currentY, sprintf('Mean diff = %s,  95%% CI [%s, %s]', ...
            smart_fmt(d), smart_fmt(ci_lo), smart_fmt(ci_hi)), 'FontSize', 10, 'Color', [0.25 0.25 0.25]);
        currentY = currentY - lineSpacing;

        % Line 2: Hedges g + percent change
        esLabel = effect_size_label(hg);
        if isfinite(pctChange)
            if pctChange < 0
                pctStr = sprintf('  |  %.1f%% reduction', abs(pctChange));
            else
                pctStr = sprintf('  |  %.1f%% increase', pctChange);
            end
        else
            pctStr = '';
        end
        if     abs(hg) >= 0.8, gCol = [0.1 0.5 0.1];
        elseif abs(hg) >= 0.5, gCol = [0.7 0.5 0];
        else,                  gCol = [0.4 0.4 0.4];
        end
        text(leftMargin + 2, currentY, sprintf('Hedges'' g = %.2f (%s effect)%s', hg, esLabel, pctStr), ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', gCol);
        currentY = currentY - lineSpacing;

        % Line 3: CI excludes zero?
        ciExcludesZero = (ci_lo > 0) || (ci_hi < 0);
        if ciExcludesZero
            text(leftMargin + 2, currentY, '95% CI excludes zero  -->  statistically significant', ...
                'FontSize', 9, 'Color', [0.1 0.55 0.2]);
        else
            text(leftMargin + 2, currentY, '95% CI includes zero  -->  not statistically significant', ...
                'FontSize', 9, 'Color', [0.5 0.5 0.5]);
        end
        currentY = currentY - lineSpacing;

        % Line 4: Direction
        if d < 0
            dirStr = sprintf('Direction: Lower in %s', expLabel);
            if isLowerBetter(metric)
                dirStr = [dirStr '  (favorable for immobilization)'];
                dirCol = [0.1 0.55 0.2];
            else
                dirCol = [0.4 0.4 0.4];
            end
        else
            dirStr = sprintf('Direction: Higher in %s', expLabel);
            if isLowerBetter(metric)
                dirStr = [dirStr '  (unfavorable)'];
                dirCol = [0.7 0.2 0.2];
            else
                dirCol = [0.4 0.4 0.4];
            end
        end
        text(leftMargin + 2, currentY, dirStr, 'FontSize', 9, 'Color', dirCol);
        currentY = currentY - lineSpacing;

        % Line 5: Precision
        ciWidth = ci_hi - ci_lo;
        if abs(d) > 1e-12
            relWidth = ciWidth / abs(d);
            if     relWidth < 1, precStr = 'Precision: HIGH (narrow CI)';     precCol = [0.1 0.5 0.1];
            elseif relWidth < 2, precStr = 'Precision: MODERATE';             precCol = [0.6 0.5 0];
            else,                precStr = 'Precision: LOW (wide CI - interpret cautiously)'; precCol = [0.7 0.3 0.1];
            end
        else
            precStr = 'Precision: Effect near zero'; precCol = [0.5 0.5 0.5];
        end
        text(leftMargin + 2, currentY, precStr, 'FontSize', 9, 'Color', precCol);
        currentY = currentY - lineSpacing;

        % Line 6: P-values
        if isfinite(pMW)
            pStr = sprintf('p-values:  Mann-Whitney %s,  Permutation %s,  FDR-corrected %s', ...
                p_fmt(pMW), p_fmt(pPerm), p_fmt(pFDR));
        else
            pStr = sprintf('p-values:  Permutation %s,  FDR-corrected %s', p_fmt(pPerm), p_fmt(pFDR));
        end
        text(leftMargin + 2, currentY, pStr, 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

    else
        text(leftMargin + 2, currentY, 'Insufficient data for analysis', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
    end

    currentY = currentY - metricGap;
end

% ==================== SUMMARY ====================
currentY = currentY + metricGap - 5;
line([leftMargin 98], [currentY currentY], 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
currentY = currentY - 3.5;

text(leftMargin, currentY, 'OVERALL SUMMARY', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.4]);
currentY = currentY - lineSpacing - 0.5;

sumCol = ternary(nSig >= nM/2, [0.1 0.55 0.2], [0.3 0.3 0.3]);
text(leftMargin + 2, currentY, ...
    sprintf('Significant effects:  %d / %d metrics  (FDR-corrected q < 0.05)', nSig, nM), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', sumCol);
currentY = currentY - lineSpacing;

nNegl = nM - nLarge - nMed - nSmall;
text(leftMargin + 2, currentY, ...
    sprintf('Effect sizes:  %d large,  %d medium,  %d small,  %d negligible', nLarge, nMed, nSmall, nNegl), ...
    'FontSize', 10, 'Color', [0.3 0.3 0.3]);
currentY = currentY - lineSpacing;

if ~isempty(pctReductions)
    text(leftMargin + 2, currentY, ...
        sprintf('Average reduction in motion metrics:  %.1f%%', mean(pctReductions)), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1]);
    currentY = currentY - lineSpacing;
end

% Footer
currentY = currentY - 2;
text(leftMargin, currentY, ...
    'Effect size benchmarks (Cohen, 1988):  |g| < 0.2 negligible,  0.2-0.5 small,  0.5-0.8 medium,  >= 0.8 large', ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45]);
currentY = currentY - lineSpacing;
text(leftMargin, currentY, ...
    sprintf('All effects computed as:  %s  minus  %s.   Negative values = less motion in experimental group.', expLabel, ctrlLabel), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45]);

hold off;

% Resize figure to fit content
usedFrac = (115 - currentY + 10) / 120;
actualH  = max(500, min(figH, figH * usedFrac * 1.15));
set(fig, 'Position', [50 50 figW actualH]);

save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% =================== RATE COMPARISON (DRIFT vs SPEED) ===================
%% ========================================================================

function plot_rate_comparison(MC, ME, pathBase, c1, c2, ctrlLabel, expLabel, dpi, vec)
mu = char(181);

if ~ismember('Drift_rate_um_per_s', MC.Properties.VariableNames) || ...
   ~ismember('Mean_speed_um_per_s', MC.Properties.VariableNames)
    fprintf('  Skipped: rate comparison (missing data)\n');
    return;
end

drC = MC.Drift_rate_um_per_s;  msC = MC.Mean_speed_um_per_s;
drE = ME.Drift_rate_um_per_s;  msE = ME.Mean_speed_um_per_s;

validC = isfinite(drC) & isfinite(msC);
validE = isfinite(drE) & isfinite(msE);

drC = drC(validC); msC = msC(validC);
drE = drE(validE); msE = msE(validE);
nC = numel(drC);   nE = numel(drE);

if nC == 0 && nE == 0
    fprintf('  Skipped: rate comparison (no valid data)\n');
    return;
end

fig = figure('Color', 'w', 'Position', [100 100 720 500], 'Visible', 'off');
hold on;

if nC > 0
    scatter(msC, drC, 80, c1, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.8);
end
if nE > 0
    scatter(msE, drE, 80, c2, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.8);
end

allMS = [msC; msE];
allDR = [drC; drE];
xMax = max(allMS) * 1.15;  if isempty(xMax) || xMax == 0, xMax = 1; end
yMax = max(allDR) * 1.15;  if isempty(yMax) || yMax == 0, yMax = 1; end

lineMax = min(xMax, yMax);
plot([0 lineMax], [0 lineMax], 'k--', 'LineWidth', 1.5);

xlim([0 xMax]);
ylim([0 yMax]);

xlabel(['Mean Speed (' mu 'm/s)'], 'FontSize', 12, 'FontWeight', 'bold');
ylabel(['Drift Rate (' mu 'm/s)'], 'FontSize', 12, 'FontWeight', 'bold');
title('Movement Pattern Analysis', 'FontSize', 14, 'FontWeight', 'bold');

h1 = scatter(nan, nan, 80, c1, 'filled', 'MarkerEdgeColor', 'w');
h2 = scatter(nan, nan, 80, c2, 'filled', 'MarkerEdgeColor', 'w');
h3 = plot(nan, nan, 'k--', 'LineWidth', 1.5);
leg = legend([h1 h2 h3], {sprintf('%s (n=%d)', ctrlLabel, nC), sprintf('%s (n=%d)', expLabel, nE), 'Identity'}, ...
    'Location', 'northeastoutside', 'FontSize', 10);
set(leg, 'Box', 'off');

annotation('textbox', [0.55 0.15 0.35 0.12], ...
    'String', {'Below line: jittering/wandering', 'On line: straight drift'}, ...
    'FontSize', 9, 'FontName', 'Arial', 'EdgeColor', [0.7 0.7 0.7], ...
    'BackgroundColor', [1 1 1 0.9], 'FitBoxToText', 'on');

set(gca, 'FontSize', 11, 'FontName', 'Arial', 'Box', 'off', 'LineWidth', 1.2, 'TickDir', 'out');

hold off;
save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ======================== RANKED SAMPLE PLOTS ===========================
%% ========================================================================

function plot_ranked_samples(MC, ME, field, metricLabel, pathBase, c1, c2, ctrlLabel, expLabel, dpi, vec)
if ~ismember(field, MC.Properties.VariableNames) || ~ismember(field, ME.Properties.VariableNames)
    fprintf('  Skipped: ranked %s (not found)\n', field);
    return;
end

allFiles  = [MC.File; ME.File];
allGroups = [MC.Group; ME.Group];
allValues = [MC.(field); ME.(field)];

validIdx  = isfinite(allValues);
allFiles  = allFiles(validIdx);
allGroups = allGroups(validIdx);
allValues = allValues(validIdx);

nSamples = numel(allValues);
if nSamples == 0
    fprintf('  Skipped: ranked %s (no data)\n', field);
    return;
end

[sortedVals, sortIdx] = sort(allValues, 'ascend');
sortedFiles  = allFiles(sortIdx);
sortedGroups = allGroups(sortIdx);
isControl    = strcmp(sortedGroups, ctrlLabel);

figHeight = min(800, max(350, 60 + nSamples * 22));
fig = figure('Color', 'w', 'Position', [100 100 820 figHeight], 'Visible', 'off');

ax = axes('Position', [0.30 0.12 0.55 0.78]);
hold(ax, 'on');

barWidth = 0.7;
for i = 1:nSamples
    if isControl(i), barColor = c1; else, barColor = c2; end
    rectangle('Position', [0, i - barWidth/2, sortedVals(i), barWidth], ...
              'FaceColor', barColor, 'EdgeColor', 'none');
end

yLabels = cell(nSamples, 1);
for i = 1:nSamples
    fname = char(sortedFiles(i));
    fname = regexprep(fname, '\.avi$', '', 'ignorecase');
    if length(fname) > 20, fname = [fname(1:17) '...']; end
    yLabels{i} = fname;
end
set(ax, 'YTick', 1:nSamples, 'YTickLabel', yLabels, 'FontSize', 7);

xlim([0, max(sortedVals) * 1.12]);
ylim([0.3, nSamples + 0.7]);
set(ax, 'YDir', 'reverse');

ctrlVals = allValues(strcmp(allGroups, ctrlLabel));
expVals  = allValues(strcmp(allGroups, expLabel));
if ~isempty(ctrlVals), xline(ax, median(ctrlVals), '--', 'Color', c1*0.7, 'LineWidth', 2); end
if ~isempty(expVals),  xline(ax, median(expVals),  '--', 'Color', c2*0.7, 'LineWidth', 2); end

nCtrl = sum(isControl);
nExp  = sum(~isControl);
h1 = patch(nan, nan, c1, 'EdgeColor', 'none');
h2 = patch(nan, nan, c2, 'EdgeColor', 'none');
leg = legend([h1 h2], {sprintf('%s (n=%d)', ctrlLabel, nCtrl), sprintf('%s (n=%d)', expLabel, nExp)}, ...
    'Location', 'northeastoutside', 'FontSize', 9);
set(leg, 'Box', 'off');

xlabel(metricLabel, 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Samples Ranked by %s', metricLabel), 'FontSize', 12, 'FontWeight', 'bold');

topHalf = floor(nSamples / 2);
if topHalf > 0
    expInTop  = sum(~isControl(1:topHalf));
    ctrlInTop = sum(isControl(1:topHalf));
    subtitle(sprintf('Lower = Better | Top half: %d %s, %d %s', ...
        ctrlInTop, ctrlLabel, expInTop, expLabel), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
end

set(ax, 'FontName', 'Arial', 'Box', 'off', 'LineWidth', 1, 'TickDir', 'out');
hold(ax, 'off');
save_figure(fig, [pathBase '.png'], dpi, vec, pathBase);
end

%% ========================================================================
%% ======================= RANKED CSV EXPORT ==============================
%% ========================================================================

function create_ranked_csv(MC, ME, csvPath, ctrlLabel, expLabel)
metrics     = {'Mean_speed_um_per_s', 'Median_step_um', 'Jitter_ratio'};
metricNames = {'MeanSpeed',           'MedianStep',     'JitterRatio'};

allFiles  = [MC.File; ME.File];
allGroups = [MC.Group; ME.Group];
nSamples  = numel(allFiles);
if nSamples == 0, return; end

T       = table();
T.File  = allFiles;
T.Group = allGroups;
allRanks = nan(nSamples, numel(metrics));

for m = 1:numel(metrics)
    if ~ismember(metrics{m}, MC.Properties.VariableNames), continue; end

    vals = [MC.(metrics{m}); ME.(metrics{m})];
    T.(metricNames{m}) = vals;
    [~, sortIdx] = sort(vals, 'ascend');
    ranks = zeros(nSamples, 1);
    ranks(sortIdx) = 1:nSamples;
    ranks(~isfinite(vals)) = NaN;
    T.([metricNames{m} '_Rank']) = ranks;
    allRanks(:, m) = ranks;
end

T.AvgRank = mean(allRanks, 2, 'omitnan');
[~, sortIdx] = sort(T.AvgRank, 'ascend');
overallRank = zeros(nSamples, 1);
overallRank(sortIdx) = 1:nSamples;
T.OverallRank = overallRank;
T = sortrows(T, 'OverallRank');

writetable(T, csvPath);
fprintf('  Saved: %s\n', csvPath);
end

%% ========================================================================
%% ============================== UTILITIES ===============================
%% ========================================================================

function save_figure(fig, pngPath, dpi, exportVec, pathBase)
drawnow;

fig.Units         = 'inches';
figPos            = fig.Position;
fig.PaperUnits    = 'inches';
fig.PaperSize     = [figPos(3) figPos(4)];
fig.PaperPosition = [0 0 figPos(3) figPos(4)];

saved = false;
try
    exportgraphics(fig, pngPath, 'Resolution', dpi, 'BackgroundColor', 'white');
    saved = true;
catch
    try
        print(fig, pngPath, '-dpng', sprintf('-r%d', dpi));
        saved = true;
    catch
    end
end

if saved
    [~, fn, ext] = fileparts(pngPath);
    fprintf('  Saved: %s%s\n', fn, ext);
else
    fprintf('  FAILED: %s\n', pngPath);
end

if exportVec && saved
    pdfPath = [pathBase '.pdf'];
    try
        exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    catch
        try
            print(fig, pdfPath, '-dpdf', '-fillpage');
        catch
        end
    end
end

close(fig);
end

function s = sanitize_filename(str)
s = lower(char(string(str)));
s = strrep(s, ' ', '_');
s = strrep(s, '+', 'plus');
s = regexprep(s, '[\\/:*?"<>|]', '-');
end

function s = smart_fmt(val)
if ~isfinite(val)
    s = 'N/A';
elseif abs(val) >= 100
    s = sprintf('%.1f', val);
elseif abs(val) >= 1
    s = sprintf('%.2f', val);
elseif abs(val) >= 0.01
    s = sprintf('%.3f', val);
elseif abs(val) >= 0.001
    s = sprintf('%.4f', val);
else
    s = sprintf('%.2g', val);
end
end

function s = p_fmt(p)
if ~isfinite(p)
    s = 'N/A';
elseif p < 0.001
    s = '<.001';
elseif p < 0.01
    s = sprintf('%.3f', p);
else
    s = sprintf('%.2f', p);
end
end

function s = effect_size_label(g)
g = abs(g);
if     g < 0.2, s = 'negligible';
elseif g < 0.5, s = 'small';
elseif g < 0.8, s = 'medium';
else,           s = 'large';
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function open_folder(folder)
try
    if ispc
        winopen(folder);
    elseif ismac
        system(['open "' folder '"']);
    else
        system(['xdg-open "' folder '" &']);
    end
catch
end
end