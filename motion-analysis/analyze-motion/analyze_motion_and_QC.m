function analyze_motion_and_QC(dataDirs, varargin)
% ANALYZE_MOTION_AND_QC  Rigid motion tracking and quality control for
%   wide-field fluorescence video recordings of Drosophila larvae.
%
%   Uses NoRMCorre rigid registration to estimate frame-by-frame
%   translational shifts, computes motion metrics, and applies automated
%   quality control flags.
%
%   USAGE:
%     analyze_motion_and_QC()                    % Opens folder picker
%     analyze_motion_and_QC('/path/to/videos')   % Process one folder
%     analyze_motion_and_QC({'/A', '/B'})         % Process multiple folders
%
%     cfg = motion_analysis_config();             % Get default config
%     cfg.pixel_size_um = 0.5;                    % Modify as needed
%     cfg.fps = 2;
%     analyze_motion_and_QC('/path', cfg)          % Use custom config
%
%   INPUTS:
%     dataDirs - (optional) Char, string, or cell array of folder paths
%                containing video files. If omitted, a folder picker opens.
%     cfg      - (optional) Configuration struct from motion_analysis_config().
%                If omitted, defaults are used. Unknown fields trigger a warning.
%
%   OUTPUTS (written to each data folder):
%     Per-file outputs:
%       *_shifts.csv             Per-frame shift estimates
%       *_motion_trace.png/pdf   Motion trace plots
%       *_rescue_meta.json       Registration rescue event log
%       *_params.json            Effective per-file registration parameters
%       *_pegged_runs.csv        Contiguous ceiling-hit episodes
%
%     Per-folder outputs:
%       motion_summary.csv       Summary statistics across files
%       motion_basic.csv         Basic motion metrics
%       motion_qc.csv            QC flags and metrics
%       motion_enriched.csv      Extended metrics (PSD, bursts, etc.)
%       motion_group_summary.csv Group-level aggregates
%       group_timecourse.csv     Per-minute motion timecourse
%       params_motion_1fps.json  Full parameter provenance record
%       thresholds.tsv           QC threshold documentation
%       column_dictionary_motion.csv  Column definitions
%       qc_overview.png/pdf      QC summary plots
%       qc_status.png/pdf        Pass/fail status strip
%       movement_overview.png/pdf  Motion overview montage
%       movement_spectrum_overview.png/pdf  Group PSD
%       motion_run.log           Processing log
%
%   CONFIGURATION:
%     All parameters (pixel size, frame rate, file pattern, QC thresholds,
%     NoRMCorre settings) are set via the config struct. See
%     motion_analysis_config() for the full list and documentation.
%
%   DEPENDENCIES:
%     - MATLAB R2020b or later (recommended: R2024b)
%     - NoRMCorre (https://github.com/flatironinstitute/NoRMCorre)
%     - Statistics and Machine Learning Toolbox (optional; for prctile)
%     - Signal Processing Toolbox (optional; for medfilt1, pwelch)
%
%   SHIFT CONVENTION (CSV output):
%     dx_px > 0 = sample moved RIGHT
%     dy_px > 0 = sample moved DOWN
%     This is MOTION direction, intuitive for visualization.
%     The pipeline applies [-dx, -dy] for correction.
%
%   See also: motion_analysis_config

% ---------- Parse config ----------
if nargin >= 2 && isstruct(varargin{1})
    cfg = merge_with_defaults(varargin{1});
elseif nargin >= 2
    error(['Unexpected second argument (type: %s). Pass a config struct from\n' ...
           'motion_analysis_config(), or call with just the data folder path.'], ...
           class(varargin{1}));
else
    cfg = motion_analysis_config();
end

% ---------- Fail-fast toolbox checks ----------
assert(exist('VideoReader','class') == 8, ...
    'VideoReader not found. A standard MATLAB installation is required.');

% ---------- Resolve folders ----------
rng(cfg.rng_seed, 'twister');

if nargin < 1 || isempty(dataDirs)
    dataDirs = prompt_for_folders();
end

% Check NoRMCorre availability early
assert(exist('NoRMCorreSetParms','file') == 2 && exist('normcorre','file') == 2, ...
    ['NoRMCorre not found on MATLAB path. Install from:\n' ...
     '  https://github.com/flatironinstitute/NoRMCorre\n' ...
     'Then add it: addpath(genpath(''/path/to/NoRMCorre''))']);

% Normalize to cellstr of char paths
if isstring(dataDirs), dataDirs = cellstr(dataDirs); end
if ~iscell(dataDirs),  dataDirs = {dataDirs};        end
dataDirs = cellfun(@(p) char(string(p)), dataDirs, 'UniformOutput', false);

% Filter non-existing with warning (continue with valid ones)
ok = true(size(dataDirs));
for i = 1:numel(dataDirs)
    if ~isfolder(dataDirs{i})
        warning('Folder does not exist and will be skipped: %s', dataDirs{i});
        ok(i) = false;
    end
end
dataDirs = dataDirs(ok);
if isempty(dataDirs)
    error('No valid folders to process.');
end

% ---------- Discovery + run ----------
N = numel(dataDirs);
for i = 1:N
    dataDir = dataDirs{i};

    % File discovery using configured pattern
    files = dir(fullfile(dataDir, cfg.file_pattern));
    fprintf('[discovery] Found %d video files in:\n  %s\n', numel(files), dataDir);
    for k = 1:numel(files)
        if k == 6 && numel(files) > 6
            fprintf('   ... (%d more)\n', numel(files)-5);
            break
        else
            fprintf('   %s\n', fullfile(files(k).folder, files(k).name));
        end
    end
    fprintf('\n=== Processing folder [%d/%d]: %s (%d files) ===\n\n', ...
        i, N, dataDir, numel(files));

    % Analyze + QC for this folder
    analyze_motion_and_QC_core(dataDir, cfg);
end

fprintf('\nAll done: analyzer + QC + enriched complete for %d folder(s).\n', N);
end

% =====================================================================
% =============== CORE ENGINE (with logging + improvements) ===========
% =====================================================================
function analyze_motion_and_QC_core(dataDir, cfg)
% One-shot: run rigid NoRMCorre + QC + enriched on video files in dataDir.

% -------- Open log file --------
logPath = fullfile(dataDir, 'motion_run.log');
logFid = fopen(logPath, 'a');
if logFid < 0
    warning('Could not open log file: %s', logPath);
    logFid = -1;
end

logf(logFid, '\n════════════════════════════════════════════════════════════════\n');
logf(logFid, 'GFP Motion Analysis Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
logf(logFid, 'Folder: %s\n', dataDir);
logf(logFid, 'Config: pixel_size_um=%.4f, fps=%g, file_pattern=%s\n', ...
    cfg.pixel_size_um, cfg.fps, cfg.file_pattern);
logf(logFid, '════════════════════════════════════════════════════════════════\n\n');

try
    % ---------- Run analyzer ----------
    [defaults_struct, versions_struct] = run_analyzer_movement(dataDir, cfg, logFid);

    % ---------- QC + Enriched ----------
    qc_and_enriched_motion_folder_movement(dataDir, cfg, defaults_struct, versions_struct, logFid);

    logf(logFid, '\n════════════════════════════════════════════════════════════════\n');
    logf(logFid, 'GFP Motion Analysis Completed: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    logf(logFid, '════════════════════════════════════════════════════════════════\n');

catch ME
    logf(logFid, '\n[ERROR] %s\n', ME.message);
    logf(logFid, 'Stack trace:\n');
    for k = 1:numel(ME.stack)
        logf(logFid, '  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    if logFid > 0, fclose(logFid); end
    rethrow(ME);
end

if logFid > 0, fclose(logFid); end
end

% =====================================================================
% ===================== ANALYZER (MOVEMENT) ============================
% =====================================================================
function [defaults_out, versions_out] = run_analyzer_movement(dataDir, cfg, logFid)

% Build defaults from config
defs = movement_defaults(cfg);
pxSize_um = defs.px_um;
fps       = defs.fps;

% Report defaults outward for JSON provenance
defaults_out = defs;

% Label style
umtxt = get_um_label(cfg.use_unicode_labels);

% Ensure NoRMCorre is on path
addpath(genpath('NoRMCorre'));
assert(exist('NoRMCorreSetParms','file') == 2 && exist('normcorre','file') == 2, ...
    'NoRMCorre not found on path. Add it first.');

% Versions (MATLAB + NoRMCorre git short hash)
versions_out = detect_versions();

logf(logFid, 'MATLAB: %s, NoRMCorre: %s\n', versions_out.MATLAB_release, versions_out.NoRMCorre_commit);

% File discovery using configured pattern
files = dir(fullfile(dataDir, cfg.file_pattern));
if isempty(files)
    warning('No files found matching %s', fullfile(dataDir, cfg.file_pattern));
    return
end

allSummary = {};
for f = 1:numel(files)
    aviPath = fullfile(files(f).folder, files(f).name);
    [~, stem] = fileparts(aviPath);
    outCSV = fullfile(files(f).folder, [stem '_shifts.csv']);

    logf(logFid, '─── %s ───\n', aviPath);

    % Guard VideoReader creation
    try
        vr = VideoReader(aviPath);
    catch ME
        logf(logFid, '[WARN] VideoReader failed for %s: %s\n', aviPath, ME.message);
        continue
    end

    % Reuse existing results if present
    if cfg.skip_if_present && exist(outCSV,'file') == 2
        logf(logFid, 'Found existing shifts CSV; skipping registration: %s\n', outCSV);
        S = compute_summary_from_csv(outCSV, pxSize_um, fps, NaN);
        if numel(S) ~= 21
            error('Summary row for %s has %d columns (expected 21).', files(f).name, numel(S));
        end
        allSummary(end+1,:) = S; %#ok<AGROW>
        continue;
    end

    H = vr.Height; W = vr.Width;

    % Per-file tuning knobs (start from config defaults)
    this_init_batch = defs.init_batch_default;
    this_bin_width  = defs.bin_width_default;
    this_us_fac     = defs.us_fac_default;

    doCrop = false; doRescale = false;
    Wuse = W; Huse = H;

    nF_est = floor(vr.Duration * max(1, vr.FrameRate));
    logf(logFid, 'Estimated frames: %d; size used: %dx%d\n', nF_est, Huse, Wuse);

    % Max shift (track initial for provenance)
    max_shift_px_initial = compute_max_shift_px(Huse, Wuse, pxSize_um, cfg.max_shift_mode, cfg.max_shift_value);
    max_shift_px = max_shift_px_initial;

    % Bounded rescue cap (for moving larvae that can traverse ~50% FOV)
    rescue_cap_px = [floor(Huse/2)-1, floor(Wuse/2)-1];
    if isscalar(max_shift_px_initial)
        rescue_cap_px = min(rescue_cap_px);
    end

    maxShift_um = max_shift_px * pxSize_um; %#ok<NASGU>
    fmt = @(v) mat2str(v);
    logf(logFid, '  Using max_shift=%s px (~%s %s) [rescue cap=%s px]\n', ...
        fmt(max_shift_px), fmt(max_shift_px.*pxSize_um), umtxt, fmt(rescue_cap_px));

    % Track rescue events
    rescue_count = 0;
    rescue_history = {};
    unstable_tracking = false;

    % Short-clip mode detection
    is_short_clip = nF_est < 300;
    if is_short_clip
        logf(logFid, '  Short clip detected (<300 frames), adjusting parameters...\n');
        this_bin_width = max(50, round(this_bin_width/2));
        this_init_batch = max(60, round(this_init_batch/2));
        this_us_fac = max(25, round(this_us_fac/2));
    end

    % RAM caps
    targetGB = cfg.target_gb;
    bytesPerFrame   = Huse * Wuse * 4; % single precision
    maxChunkByRAM   = max(120, floor((targetGB*1e9) / bytesPerFrame));
    chunk_frames    = min(defs.chunk_frames_default, maxChunkByRAM);
    if this_init_batch > chunk_frames
        logf(logFid, '  Capping init_batch from %d to %d for RAM budget (~%.1f GB target)\n', this_init_batch, chunk_frames, targetGB);
        this_init_batch = chunk_frames;
    end
    if chunk_frames < defs.chunk_frames_default
        logf(logFid, '  Capping chunk_frames from %d to %d for RAM budget (~%.1f GB target)\n', defs.chunk_frames_default, chunk_frames, targetGB);
    end

    % NoRMCorre options
    options_mc = NoRMCorreSetParms('d1',Huse,'d2',Wuse, ...
        'bin_width',this_bin_width, ...
        'max_shift',max_shift_px, ...
        'us_fac',this_us_fac, ...
        'init_batch',this_init_batch);

    % Initial template
    logf(logFid, 'Building initial template from first %d frames (~%.1f min) ...\n', this_init_batch, this_init_batch/fps/60);
    nInit = 0;
    initBuf = zeros(Huse, Wuse, this_init_batch, 'single');
    while hasFrame(vr) && nInit < this_init_batch
        nInit = nInit + 1;
        fr = readFrame_gray(vr);
        if doCrop, error('Crop map not set in this run.'); end %#ok<UNRCH>
        fr = single(fr);
        if doRescale, fr = rescale01(fr); end
        initBuf(:,:,nInit) = fr;
    end
    initBuf = initBuf(:,:,1:nInit);
    if nInit == 0, error('No frames read. Is the video file valid/uncompressed?'); end

    template = []; shifts_init = [];
    try
        [~, shifts_init, template] = normcorre(initBuf, options_mc); %#ok<ASGLU>
    catch ME
        if contains(ME.message,'Out of memory','IgnoreCase',true)
            logf(logFid, '  OOM during template; retry with lighter settings (us_fac/2, init_batch/2)...\n');
            this_us_fac = max(10, round(this_us_fac/2));
            this_init_batch = max(120, round(this_init_batch/2));
            options_mc = NoRMCorreSetParms('d1',Huse,'d2',Wuse, ...
                'bin_width',this_bin_width,'max_shift',max_shift_px, ...
                'us_fac',this_us_fac,'init_batch',this_init_batch);
            initBuf = initBuf(:,:,1:min(size(initBuf,3), this_init_batch));
            [~, shifts_init, template] = normcorre(initBuf, options_mc); %#ok<ASGLU>
        else
            rethrow(ME);
        end
    end
    clear initBuf;
    [dx_all, dy_all] = extract_dx_dy(shifts_init);

    % Streaming registration
    options_stream = options_mc; options_stream.init_batch = 0;
    chunk = zeros(Huse, Wuse, chunk_frames, 'single'); nIn = 0;
    chunk_rescue_count = 0;

    while hasFrame(vr)
        fr = readFrame_gray(vr);
        if doRescale, fr = rescale01(fr); end
        nIn = nIn + 1; chunk(:,:,nIn) = fr;

        if nIn == chunk_frames || ~hasFrame(vr)
            Yc = chunk(:,:,1:nIn);
            chunk_rescue_count = 0;

            % OOM ladder
            attempts = 0;
            while true
                try
                    [~, shifts_c, template] = normcorre(Yc, options_stream, template); %#ok<ASGLU>
                    break
                catch ME
                    if ~contains(ME.message,'Out of memory','IgnoreCase',true), rethrow(ME); end
                    attempts = attempts + 1;
                    if attempts == 1 && size(Yc,3) > 180
                        newN = max(120, floor(size(Yc,3)/2));
                        logf(logFid, '  [OOM] retry with smaller chunk (%d -> %d frames)\n', size(Yc,3), newN);
                        Yc = Yc(:,:,1:newN);
                    elseif attempts == 2 && isfield(options_stream,'us_fac') && options_stream.us_fac > 10
                        options_stream.us_fac = max(10, round(options_stream.us_fac/2));
                        logf(logFid, '  [OOM] retry with lower us_fac=%d\n', options_stream.us_fac);
                    else
                        rethrow(ME);
                    end
                end
            end

            % Bounded near-ceiling rescue (up to 3 times per chunk)
            [dx_try, dy_try] = extract_dx_dy(shifts_c);

            while isfield(options_stream,'max_shift') && ~isempty(dx_try) && chunk_rescue_count < 3
                ms = options_stream.max_shift;
                if isscalar(ms), msy = ms; msx = ms; else, msy = ms(1); msx = ms(2); end
                if isscalar(rescue_cap_px), capy = rescue_cap_px; capx = rescue_cap_px; else, capy = rescue_cap_px(1); capx = rescue_cap_px(2); end

                hit = any(abs(dx_try) > 0.90*msx) || any(abs(dy_try) > 0.90*msy);
                if ~hit, break; end

                bump_y = min(capy, max(msy+10, ceil(1.5*msy)));
                bump_x = min(capx, max(msx+10, ceil(1.5*msx)));
                new_ms_vec = [bump_y, bump_x];

                if all([bump_y <= msy, bump_x <= msx]) || any([bump_y > capy, bump_x > capx])
                    break; % no meaningful growth left
                end

                try
                    opts_bump = options_stream;
                    if isscalar(ms)
                        opts_bump.max_shift = max(new_ms_vec);
                    else
                        opts_bump.max_shift = new_ms_vec;
                    end
                    logf(logFid, '  [RESCUE] near ceiling %s -> max_shift=%s px\n', ...
                        fmt(ms), fmt(opts_bump.max_shift));
                    [~, shifts_c2, template] = normcorre(Yc, opts_bump, template); %#ok<ASGLU>
                    shifts_c = shifts_c2;
                    options_stream.max_shift = opts_bump.max_shift;

                    rescue_count = rescue_count + 1;
                    chunk_rescue_count = chunk_rescue_count + 1;
                    frame_start = numel(dx_all) + 1;
                    frame_end = frame_start + size(Yc,3) - 1;
                    rescue_history{end+1} = sprintf('%s->%s px, frames %d-%d', ...
                        fmt(ms), fmt(opts_bump.max_shift), frame_start, frame_end); %#ok<AGROW>
                    [dx_try, dy_try] = extract_dx_dy(shifts_c);
                    continue
                catch
                    break
                end
            end

            % If still hitting ceiling after max rescues, mark unstable
            if chunk_rescue_count >= 3
                ms = options_stream.max_shift;
                if isscalar(ms), msy = ms; msx = ms; else, msy = ms(1); msx = ms(2); end
                if any(abs(dx_try) > 0.90*msx) || any(abs(dy_try) > 0.90*msy)
                    logf(logFid, '  [RESCUE] Hit ceiling after 3 bumps, marking UNSTABLE\n');
                    unstable_tracking = true;
                end
            end

            [dx_c, dy_c] = extract_dx_dy(shifts_c);
            dx_all = [dx_all; dx_c(:)]; %#ok<AGROW>
            dy_all = [dy_all; dy_c(:)]; %#ok<AGROW>
            nIn = 0;
        end
    end

    nF = numel(dx_all);

    % Early sanity check
    if abs(nF - nF_est) / max(1, nF_est) > 0.1
        logf(logFid, '[WARN] Frame count differs from duration estimate by >10%% (actual=%d, est=%d)\n', nF, nF_est);
    end

    % Determine the final/effective max_shift actually used
    if isfield(options_stream,'max_shift') && ~isempty(options_stream.max_shift)
        eff_ms = options_stream.max_shift;
    else
        eff_ms = max_shift_px;
    end
    if isscalar(eff_ms)
        eff_max_shift_px = eff_ms;
    else
        eff_max_shift_px = max(eff_ms);
    end
    eff_maxShift_um = eff_max_shift_px * pxSize_um;

    logf(logFid, 'Total frames processed: %d (us_fac=%d, init_batch=%d, chunk<=%d)\n', ...
        nF, this_us_fac, this_init_batch, chunk_frames);

    if unstable_tracking, unstable_tag = ' [UNSTABLE]'; else, unstable_tag = ''; end
    logf(logFid, '  Max shift: initial=%s px, effective=%s px, rescues=%d%s\n\n', ...
        fmt(max_shift_px_initial), fmt(eff_ms), rescue_count, unstable_tag);

    % Per-frame CSV
    dx_um = dx_all * pxSize_um; dy_um = dy_all * pxSize_um;
    steps_um = hypot([0; diff(dx_um)], [0; diff(dy_um)]);

    % Check if step distribution hits search window
    p99_step = Pctl(steps_um, 99);
    if p99_step > 0.9 * eff_maxShift_um
        logf(logFid, '[WARN] Step distribution (P99=%.2f %s) approaches search window (%.2f %s)\n', ...
            p99_step, umtxt, eff_maxShift_um, umtxt);
    end

    t_s   = (0:nF-1).'/fps;
    t_min = t_s/60;
    T = table((1:nF)', t_s, t_min, dx_all, dy_all, dx_um, dy_um, steps_um, ...
        repmat(max_shift_px_initial, nF,1), repmat(max_shift_px_initial*pxSize_um, nF,1), ...
        repmat(eff_max_shift_px, nF,1), repmat(eff_maxShift_um, nF,1), ...
        repmat(rescue_count, nF,1), repmat(rescue_cap_px, nF,1), ...
        repmat(logical(unstable_tracking), nF,1), ...
        repmat(fps, nF,1), repmat(pxSize_um, nF,1), ...
        'VariableNames', {'Frame','time_s','time_min','dx_px','dy_px','dx_um','dy_um','displacement_um', ...
        'max_shift_px_initial','max_shift_um_initial','max_shift_px','max_shift_um', ...
        'rescue_count','rescue_cap_px','unstable_tracking','fps_used','px_um'});

    T.MATLAB_release   = repmat(string(versions_out.MATLAB_release), nF, 1);
    T.NoRMCorre_commit = repmat(string(versions_out.NoRMCorre_commit), nF, 1);

    % Add optional [y x] columns if using vector max_shift
    if exist('eff_ms','var') && ~isscalar(eff_ms)
        T.max_shift_px_y = repmat(eff_ms(1), nF, 1);
        T.max_shift_px_x = repmat(eff_ms(2), nF, 1);
        T.max_shift_um_y = repmat(eff_ms(1) * pxSize_um, nF, 1);
        T.max_shift_um_x = repmat(eff_ms(2) * pxSize_um, nF, 1);
    end

    % Add rescue_cap columns if vector
    if ~isscalar(rescue_cap_px)
        if numel(rescue_cap_px) ~= 2
            error('rescue_cap_px must be scalar or 2-element vector, got %d elements', numel(rescue_cap_px));
        end
        T.rescue_cap_px_y = repmat(rescue_cap_px(1), nF, 1);
        T.rescue_cap_px_x = repmat(rescue_cap_px(2), nF, 1);
    end

    outCSV = fullfile(files(f).folder, [stem '_shifts.csv']);
    writetable(T, outCSV);

    % Write rescue metadata JSON
    meta = struct('file', files(f).name, ...
                  'rescue_count', rescue_count, ...
                  'rescue_cap_px', rescue_cap_px, ...
                  'unstable_tracking', logical(unstable_tracking), ...
                  'rescue_history', {rescue_history});
    metaJSON = fullfile(files(f).folder, [stem '_rescue_meta.json']);
    fid = fopen(metaJSON, 'w');
    if fid > 0
        fwrite(fid, jsonencode(meta), 'char');
        fclose(fid);
    end

    % Write per-file effective parameters JSON
    file_params = struct('bin_width', this_bin_width, ...
                        'init_batch', this_init_batch, ...
                        'us_fac', this_us_fac, ...
                        'max_shift_initial', max_shift_px_initial, ...
                        'max_shift_effective', eff_ms, ...
                        'rescue_count', rescue_count, ...
                        'unstable', logical(unstable_tracking));
    paramsJSON = fullfile(files(f).folder, [stem '_params.json']);
    fid = fopen(paramsJSON, 'w');
    if fid > 0
        fwrite(fid, jsonencode(file_params), 'char');
        fclose(fid);
    end

    % Write pegged runs table
    if exist('eff_ms','var') && ~isscalar(eff_ms)
        thr_y_um = eff_ms(1)*pxSize_um;
        thr_x_um = eff_ms(2)*pxSize_um;
    else
        thr_y_um = eff_maxShift_um;
        thr_x_um = eff_maxShift_um;
    end
    is_pegged = (abs(dx_um) > 0.98*thr_x_um) | (abs(dy_um) > 0.98*thr_y_um);
    [peg_starts, peg_ends] = find_runs(is_pegged);
    if ~isempty(peg_starts)
        peggedTable = table(peg_starts, peg_ends, (peg_ends-peg_starts+1)/fps, ...
            'VariableNames', {'start_idx','end_idx','duration_s'});
        writetable(peggedTable, fullfile(files(f).folder, [stem '_pegged_runs.csv']));
    end

    % Quick plot (enhanced with shaded regions) - PNG + PDF
    try
        set_figure_defaults();

        lim_dx   = Pctl(abs(dx_um), 99);
        lim_dy   = Pctl(abs(dy_um), 99);
        lim_step = Pctl(steps_um, 99);
        pad = 1.15;
        ydx   = max(lim_dx*pad, 5);
        ydy   = max(lim_dy*pad, 5);
        ystep = max(lim_step*pad, 0.5);

        % Spike threshold for visualization
        SPIKE_THR_RATE_UM_S = cfg.qc_spike_thr_um_s;
        spike_thr_um = SPIKE_THR_RATE_UM_S / fps;
        is_spike = steps_um > spike_thr_um;

        fign = figure('Visible','off','Color','w','Position',[200 200 980 760]);
        tEnd = t_min(end);

        subplot(3,1,1); plot(t_min, dx_um, 'LineWidth', 1); grid on;
        ylabel(['Delta x (' umtxt ')']); title(strrep(files(f).name,'_','\_'));
        ylim([-ydx ydx]); xlim([0, tEnd]);

        subplot(3,1,2); plot(t_min, dy_um, 'LineWidth', 1); grid on;
        ylabel(['Delta y (' umtxt ')']); ylim([-ydy ydy]); xlim([0, tEnd]);

        subplot(3,1,3); plot(t_min, steps_um, 'LineWidth', 1); grid on; hold on;
        shade_segments(t_min, is_pegged, [1 0.85 0.85]);
        shade_segments(t_min, is_spike,  [1 0.95 0.85]);
        ylabel(['step (' umtxt ')']); xlabel('time (min)'); ylim([0 ystep]); xlim([0, tEnd]);

        outPNG = fullfile(files(f).folder, sprintf('%s_motion_trace.png', stem));
        outPDF = fullfile(files(f).folder, sprintf('%s_motion_trace.pdf', stem));
        exportgraphics(fign, outPNG, 'Resolution', 300);
        exportgraphics(fign, outPDF);
        close(fign);
    catch
    end

    % Summary row
    win  = max(3, 2*floor((60*fps)/2)+1);
    try
        dx_s = medfilt1(dx_um, win, 'omitnan', 'truncate');
        dy_s = medfilt1(dy_um, win, 'omitnan', 'truncate');
    catch
        dx_s = movmedian(dx_um, win, 'omitnan');
        dy_s = movmedian(dy_um, win, 'omitnan');
    end
    net_drift  = hypot(dx_s(end)-dx_s(1), dy_s(end)-dy_s(1));
    duration_s = (nF - 1) / max(1, fps);
    drift_rate = net_drift / max(1e-9, duration_s);

    med_step  = median(steps_um, 'omitnan');
    p95_step  = Pctl(steps_um,95);
    cum_path  = sum(steps_um,'omitnan');

    % For GFP: whole movie stats (no settle)
    settle_frames = 0; %#ok<NASGU>
    med_step_ps   = med_step;
    p95_step_ps   = p95_step;
    cum_path_ps   = cum_path;
    net_drift_ps  = net_drift; %#ok<NASGU>
    drift_rate_ps = drift_rate; %#ok<NASGU>

    settle_min_report = 0;

    summaryRow = { ...
        files(f).name, nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, NaN, NaN, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, NaN, NaN, ...
        fps, pxSize_um, max_shift_px_initial * pxSize_um, eff_maxShift_um, settle_min_report};

    if numel(summaryRow) ~= 21
        error('Summary row for %s has %d columns (expected 21).', files(f).name, numel(summaryRow));
    end
    allSummary(end+1,:) = summaryRow; %#ok<AGROW>
end

if ~isempty(allSummary)
    expectedCols = 21;
    if size(allSummary, 2) ~= expectedCols
        error('Expected %d columns in summary but got %d. Check compute_summary_from_csv.', ...
              expectedCols, size(allSummary, 2));
    end

    summaryTable = cell2table(allSummary, ...
        'VariableNames', {'File','n_frames', ...
        'MedianDisp_um','Pct95Disp_um','CumPath_um','NetDrift_um','DriftRate_um_per_s','Pegged_frac','Flatline_flag', ...
        'MedianDisp_postSettle_um','Pct95Disp_postSettle_um','CumPath_postSettle_um','NetDrift_postSettle_um','DriftRate_postSettle_um_per_s','Pegged_postSettle_frac','Flatline_postSettle_flag', ...
        'fps_used','px_um','max_shift_um_initial','max_shift_um_effective','settle_min'});

    summaryTable.MATLAB_release   = repmat(string(versions_out.MATLAB_release), height(summaryTable), 1);
    summaryTable.NoRMCorre_commit = repmat(string(versions_out.NoRMCorre_commit), height(summaryTable), 1);

    writetable(summaryTable, fullfile(dataDir,'motion_summary.csv'));
    logf(logFid, 'Wrote summary: %s\n', fullfile(dataDir,'motion_summary.csv'));
end
end

% =====================================================================
% ================= QC + ENRICHED (GFP-friendly) ======================
% =====================================================================
function qc_and_enriched_motion_folder_movement(dataDir, cfg, defs_in, versions_in, logFid)
% Movement QC with configurable parameters. GFP-friendly: only tracking
% integrity matters for the PASS flag.

px_um_default = cfg.pixel_size_um;
fps_default   = cfg.fps;

% ---- HARD thresholds (tracking integrity ONLY for GFP) ----
max_pegged_frac   = cfg.qc_max_pegged_frac;
peg_run_min_s     = cfg.qc_peg_run_min_s;
nan_frac_max      = cfg.qc_nan_frac_max;

% ---- SOFT thresholds (descriptive only, not gating for GFP) ----
good_post_med_um   = cfg.qc_good_post_med_um;
spike_thr_rate_um_s= cfg.qc_spike_thr_um_s;
max_spike_frac     = cfg.qc_max_spike_frac;
flat_eps_um        = max(0.05, 0.25*px_um_default);

% For reporting only (not gating):
imp_jump_rate_um_s = cfg.qc_imp_jump_um_s;
imp_jump_frac_max  = cfg.qc_imp_jump_frac_max; %#ok<NASGU>

% ---- Settle policy (analyze whole movie for GFP) ----
base_settle_min    = 0;
use_auto_settle_for_QC = false; %#ok<NASGU>

umtxt = get_um_label(cfg.use_unicode_labels);

% Derive CSV discovery pattern from video file pattern
[~, pat_stem, ~] = fileparts(cfg.file_pattern);
csv_pattern = [pat_stem '_shifts.csv'];
csvs = dir(fullfile(dataDir, csv_pattern));
if isempty(csvs)
    logf(logFid, '[WARN] QC: no per-movie CSVs found matching %s in %s. Nothing to do.\n', csv_pattern, dataDir);
    return
end
if ~cfg.include_cropped
    csvs = csvs(~contains({csvs.name},'cropped','IgnoreCase',true));
end

logf(logFid, '\n─── QC PROCESSING ───\n');
logf(logFid, 'GFP mode: QC gates on tracking integrity only\n\n');

rows_basic = {}; rows_qc = {}; rows_enr = {};
tc_all = {}; psd_all = {};

for i = 1:numel(csvs)
    inCSV = fullfile(csvs(i).folder, csvs(i).name);
    [~,stem] = fileparts(inCSV);
    stem = regexprep(stem, '_shifts$', '');
    T = readtable(inCSV);

    px_um = px_um_default; if ismember('px_um', T.Properties.VariableNames), px_um = T.px_um(1); end
    fps   = fps_default;   if ismember('fps_used', T.Properties.VariableNames), fps = T.fps_used(1); end

    % Adaptive PSD bands based on Nyquist
    ny = fps/2;
    band_low  = [0.01, min(0.05, 0.2*ny)];
    band_mid  = [min(0.05, 0.2*ny), min(0.2, 0.8*ny)];
    band_high = [min(0.2, 0.8*ny), ny];

    % Per-frame thresholds derived from fps
    spike_thr_um         = spike_thr_rate_um_s / fps;  %#ok<NASGU>
    spike_min_s          = 2.0;
    spike_min_consec     = max(1, ceil(spike_min_s * fps));

    imp_jump_um          = imp_jump_rate_um_s / fps; %#ok<NASGU>
    imp_min_consec       = max(1, ceil(0.5 * fps));

    peg_run_len_fail     = max(1, ceil(peg_run_min_s * fps));

    % Try to get image dimensions from the video file
    [~, video_ext] = fileparts(cfg.file_pattern);
    aviPath = fullfile(csvs(i).folder, [stem video_ext]);
    if exist(aviPath,'file') == 2
        vr = VideoReader(aviPath); H = vr.Height; W = vr.Width;
    else
        H = 1024; W = 1024;
    end

    % Prefer the CSV's max_shift if present, else compute
    if ismember('max_shift_um', T.Properties.VariableNames)
        maxShift_um = T.max_shift_um(1);
    else
        max_shift_px = compute_max_shift_px(H, W, px_um, cfg.max_shift_mode, cfg.max_shift_value);
        maxShift_um  = max_shift_px * px_um;
    end

    % Read rescue metadata if available
    if ismember('rescue_count', T.Properties.VariableNames)
        rescue_count = T.rescue_count(1);
    else
        rescue_count = 0;
    end
    if ismember('unstable_tracking', T.Properties.VariableNames)
        unstable_tracking = logical(T.unstable_tracking(1));
    else
        unstable_tracking = false;
    end

    if all(ismember({'dx_um','dy_um'}, T.Properties.VariableNames))
        dx_um = T.dx_um; dy_um = T.dy_um;
    else
        dx_um = T.dx_px * px_um; dy_um = T.dy_px * px_um;
    end
    nF = numel(dx_um);
    steps_um = hypot([0; diff(dx_um)], [0; diff(dy_um)]);
    duration_s = (nF - 1) / max(1, fps);
    duration_min = duration_s / 60;
    drift_speed_um_s = NaN;

    % --- For GFP: analyze entire movie ---
    settle_start_min = 0;
    pre_idx  = 1:0;          % empty
    post_idx = 1:nF;         % whole movie
    tail_length_min = duration_min;

    % Still compute auto-settle for reporting only
    auto_settle_min = NaN;
    try
        sm = movmedian(steps_um, round(60*fps), 'omitnan');
        below = sm < good_post_med_um;
        runlen = 0; startIdx = NaN;
        for k=1:numel(below)
            if below(k)
                if isnan(startIdx), startIdx = k; end
                runlen = runlen + 1;
                if runlen >= round(3*60*fps)
                    auto_settle_min = (startIdx-1)/(60*fps);
                    break
                end
            else
                runlen = 0; startIdx = NaN;
            end
        end
    catch
    end

    % --- Reliability (hard) metrics over FULL movie ---
    if all(ismember({'max_shift_um_y','max_shift_um_x'}, T.Properties.VariableNames))
        thr_y_um = T.max_shift_um_y(1);
        thr_x_um = T.max_shift_um_x(1);
    else
        thr_y_um = maxShift_um;
        thr_x_um = maxShift_um;
    end
    is_pegged = (abs(dx_um) > 0.98*thr_x_um) | (abs(dy_um) > 0.98*thr_y_um);
    pegged    = mean(is_pegged);
    nan_frac  = mean(~isfinite(dx_um) | ~isfinite(dy_um));
    max_pegrun= max_run(is_pegged);

    flatline = (std(dx_um) < flat_eps_um) && (std(dy_um) < flat_eps_um);

    win  = max(3, 2*floor((60*fps)/2)+1);
    try
        dx_s = medfilt1(dx_um, win, 'omitnan', 'truncate');
        dy_s = medfilt1(dy_um, win, 'omitnan', 'truncate');
    catch
        dx_s = movmedian(dx_um, win, 'omitnan');
        dy_s = movmedian(dy_um, win, 'omitnan');
    end
    net_drift  = hypot(dx_s(end)-dx_s(1), dy_s(end)-dy_s(1));
    drift_rate = net_drift / max(1e-9, duration_s);
    drift_speed_um_s = drift_rate;

    med_step = median(steps_um,'omitnan');
    mean_step = mean(steps_um,'omitnan');
    sd_step   = std(steps_um, 'omitnan');
    p95_step  = Pctl(steps_um,95);
    p99_step  = Pctl(steps_um,99);
    cum_path  = sum(steps_um,'omitnan');

    drift_dx = dx_s(end)-dx_s(1);
    drift_dy = dy_s(end)-dy_s(1);
    drift_angle_deg = atan2d(drift_dy, drift_dx);

    pct_below_0p10 = mean(steps_um <= 0.10, 'omitnan');
    pct_below_0p20 = mean(steps_um <= 0.20, 'omitnan');

    med_step_pre = NaN;  % No pre-settle for GFP

    % --- Post-tail arrays (whole movie for GFP) ---
    if numel(post_idx) >= 2
        steps_post   = steps_um(post_idx);
        step_rate_post = steps_post * fps;

        imp_jump_frac_ps = frac_above_threshold(step_rate_post, imp_jump_rate_um_s, imp_min_consec);
        spike_frac_ps    = frac_above_threshold(step_rate_post, spike_thr_rate_um_s, spike_min_consec);

        med_step_ps  = median(steps_post,'omitnan');
        mean_step_ps = mean(steps_post,'omitnan');
        p95_step_ps  = Pctl(steps_post,95);
        p99_step_ps  = Pctl(steps_post,99);
        cum_path_ps  = sum(steps_post,'omitnan');

        net_drift_ps  = net_drift;
        drift_rate_ps = drift_rate;

        pegged_ps   = pegged;
        flatline_ps = flatline;

        jitter_ratio = cum_path_ps / max(1e-9, net_drift_ps);
        jitter_ratio_display = min(jitter_ratio, 1000);

        [f_psd, P] = compute_psd_safe(steps_post, fps);
        [dom_f, dom_pow]   = dominant_peak(f_psd, P, [0.01 0.50]);
        dom_period_min     = (dom_f>0) * (1/dom_f) / 60;
        frac_low  = band_power_frac(f_psd, P, band_low);
        frac_mid  = band_power_frac(f_psd, P, band_mid);
        frac_high = band_power_frac(f_psd, P, band_high);
        psd_all{end+1} = struct('f',f_psd,'P',P); %#ok<AGROW>
    else
        [med_step_ps,mean_step_ps,p95_step_ps,p99_step_ps,cum_path_ps,net_drift_ps,drift_rate_ps,pegged_ps,flatline_ps,spike_frac_ps,imp_jump_frac_ps,jitter_ratio,jitter_ratio_display] = deal(NaN);
        [dom_f, dom_pow, dom_period_min, frac_low, frac_mid, frac_high] = deal(NaN);
    end

    settle_ratio = NaN;  % Not meaningful for GFP

    % ---- HARD PASS: tracking integrity ONLY (GFP) ----
    failReasons = strings(0,1);
    if pegged   >= max_pegged_frac,   failReasons(end+1) = "pegged_frac";      end
    if max_pegrun >= peg_run_len_fail,failReasons(end+1) = "ceiling_run";      end
    if nan_frac >= nan_frac_max,      failReasons(end+1) = "nan_frac";         end
    if flatline,                      failReasons(end+1) = "flatline";         end
    if unstable_tracking,             failReasons(end+1) = "unstable_tracking"; end

    qc_pass_tracking = isempty(failReasons);
    qc_soft_still    = (isfinite(med_step_ps) && med_step_ps < good_post_med_um) ...
                       && (isfinite(spike_frac_ps) && spike_frac_ps < max_spike_frac);

    qc_pass = qc_pass_tracking;

    if isempty(failReasons), qc_reason = "OK";
    else, qc_reason = strjoin(failReasons, ",");
    end

    HealthFlag = "OK";
    if unstable_tracking
        HealthFlag = "UNSTABLE_TRACKING";
    elseif flatline
        HealthFlag = "FLATLINE";
    elseif ~qc_pass_tracking
        if nan_frac >= nan_frac_max, HealthFlag = "NAN_EXCESS";
        elseif max_pegrun >= peg_run_len_fail, HealthFlag = "CEILING_RUN";
        elseif pegged >= max_pegged_frac, HealthFlag = "PEGGED";
        end
    end

    % -------- BASIC --------
    rows_basic(end+1,:) = { ...
        [stem video_ext], nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, pegged, flatline, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, pegged_ps, flatline_ps, ...
        fps, px_um, maxShift_um, settle_start_min }; %#ok<AGROW>

    % -------- QC --------
    rows_qc(end+1,:) = { ...
        [stem video_ext], nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, pegged, flatline, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, pegged_ps, flatline_ps, ...
        fps, px_um, maxShift_um, settle_start_min, ...
        med_step_pre, settle_ratio, spike_frac_ps, jitter_ratio_display, ...
        logical(qc_pass), ...
        tail_length_min, auto_settle_min, ...
        logical(qc_pass_tracking), logical(qc_soft_still), string(qc_reason), string(HealthFlag), ...
        nan_frac, rescue_count, logical(unstable_tracking)}; %#ok<AGROW>

    % -------- ENRICHED --------
    rows_enr(end+1,:) = { ...
        [stem video_ext], nF, ...
        mean_step, sd_step, p99_step, ...
        mean_step_ps, p99_step_ps, ...
        pct_below_0p10, pct_below_0p20, ...
        drift_angle_deg, auto_settle_min, ...
        burst_count(steps_post_if(post_idx,steps_um), 0.2, fps), burst_rate(steps_post_if(post_idx,steps_um), 0.2, fps), ...
        burst_count(steps_post_if(post_idx,steps_um), 0.5, fps), burst_rate(steps_post_if(post_idx,steps_um), 0.5, fps), ...
        burst_count(steps_post_if(post_idx,steps_um), 1.0, fps), burst_rate(steps_post_if(post_idx,steps_um), 1.0, fps), ...
        dom_f, dom_period_min, dom_pow, ...
        frac_low, frac_mid, frac_high, ...
        med_step, med_step_ps, cum_path, cum_path_ps, ...
        pegged, pegged_ps, spike_frac_ps, jitter_ratio_display, drift_speed_um_s, ...
        logical(qc_pass), ...
        settle_start_min, tail_length_min, ...
        logical(qc_pass_tracking), logical(qc_soft_still), string(HealthFlag), ...
        nan_frac, rescue_count, logical(unstable_tracking)}; %#ok<AGROW>

    [min_idx, per_min_median] = per_minute_median(steps_um, fps);
    tc_all{end+1} = struct('minute',min_idx,'file_median',per_min_median); %#ok<AGROW>
end

% Headers
hdr_basic = {'File','n_frames', ...
    'MedianDisp_um','Pct95Disp_um','CumPath_um','NetDrift_um','DriftRate_um_per_s','Pegged_frac','Flatline_flag', ...
    'MedianDisp_postSettle_um','Pct95Disp_postSettle_um','CumPath_postSettle_um','NetDrift_postSettle_um','DriftRate_postSettle_um_per_s','Pegged_postSettle_frac','Flatline_postSettle_flag', ...
    'fps_used','px_um','max_shift_um','settle_min'};

hdr_qc = [hdr_basic, ...
    {'MedianDisp_preSettle_um','SettleRatio_post_over_pre','SpikeFrac_postSettle_over_0p5um_per_s','JitterRatio_post_cumPath_over_netDrift','QC_PASS', ...
    'TailLength_min','AutoSettle_min', ...
    'QC_PASS_tracking','QC_soft_still','QC_reason','HealthFlag','nan_frac','rescue_count','unstable_tracking'}];

hdr_enr = {'File','n_frames', ...
    'MeanStep_um','SDStep_um','P99Step_um', ...
    'MeanStep_postSettle_um','P99Step_postSettle_um', ...
    'PctTime_step_le_0p10um','PctTime_step_le_0p20um', ...
    'DriftAngle_deg','AutoSettle_min', ...
    'Bursts_post_0p2um','Rate_post_0p2um_per_min', ...
    'Bursts_post_0p5um','Rate_post_0p5um_per_min', ...
    'Bursts_post_1p0um','Rate_post_1p0um_per_min', ...
    'PSD_DominantFreq_Hz','PSD_DominantPeriod_min','PSD_DominantPower', ...
    'PSD_FracLow','PSD_FracMid','PSD_FracHigh', ...
    'MedianStep_um','MedianStep_postSettle_um', ...
    'CumPath_um','CumPath_postSettle_um', ...
    'Pegged_frac','Pegged_postSettle_frac', ...
    'SpikeFrac_postSettle_over_0p5um_per_s','JitterRatio_post_cumPath_over_netDrift','DriftSpeed_um_s', ...
    'QC_PASS', ...
    'TailStart_min','TailLength_min', ...
    'QC_PASS_tracking','QC_soft_still','HealthFlag','nan_frac','rescue_count','unstable_tracking'};

perfile_basic = cell2table(rows_basic,'VariableNames',hdr_basic);
perfile_qc    = cell2table(rows_qc,   'VariableNames',hdr_qc);
perfile_enr   = cell2table(rows_enr,  'VariableNames',hdr_enr);

% Append versions
if ~isempty(fieldnames(versions_in))
    perfile_basic.MATLAB_release   = repmat(string(versions_in.MATLAB_release), height(perfile_basic), 1);
    perfile_basic.NoRMCorre_commit = repmat(string(versions_in.NoRMCorre_commit), height(perfile_basic), 1);

    perfile_qc.MATLAB_release   = repmat(string(versions_in.MATLAB_release), height(perfile_qc), 1);
    perfile_qc.NoRMCorre_commit = repmat(string(versions_in.NoRMCorre_commit), height(perfile_qc), 1);

    perfile_enr.MATLAB_release   = repmat(string(versions_in.MATLAB_release), height(perfile_enr), 1);
    perfile_enr.NoRMCorre_commit = repmat(string(versions_in.NoRMCorre_commit), height(perfile_enr), 1);
end

out_basic   = fullfile(dataDir,'motion_basic.csv');
out_qc      = fullfile(dataDir,'motion_qc.csv');
out_enr     = fullfile(dataDir,'motion_enriched.csv');
writetable(perfile_basic, out_basic);
writetable(perfile_qc,    out_qc);
writetable(perfile_enr,   out_enr);
logf(logFid, 'QC/ENR: wrote %s, %s, %s\n', out_basic, out_qc, out_enr);

% Group aggregates
G = table;
G.n_files        = height(perfile_qc);
G.n_pass         = sum(perfile_qc.QC_PASS);
G.pass_rate      = G.n_pass / max(1,G.n_files);
G.n_pass_tracking = sum(perfile_qc.QC_PASS_tracking);
G.n_soft_still   = sum(perfile_qc.QC_soft_still);
m = @(x) mean(x,'omitnan'); s = @(x) std(x,'omitnan');
G.mean_MedianDisp_post_um   = m(perfile_qc.MedianDisp_postSettle_um);
G.sd_MedianDisp_post_um     = s(perfile_qc.MedianDisp_postSettle_um);
G.mean_SettleRatio          = m(perfile_qc.SettleRatio_post_over_pre);
G.sd_SettleRatio            = s(perfile_qc.SettleRatio_post_over_pre);
G.mean_AutoSettle_min       = m(perfile_qc.AutoSettle_min);
G.mean_TailLength_min       = m(perfile_qc.TailLength_min);
writetable(G, fullfile(dataDir,'motion_group_summary.csv'));
logf(logFid, 'QC/ENR: wrote motion_group_summary.csv\n');

% Log QC summary
logf(logFid, '\n═══ QC SUMMARY ═══\n');
logf(logFid, 'Total files: %d\n', G.n_files);
logf(logFid, '  QC_PASS (tracking): %d (%.1f%%)\n', G.n_pass, G.pass_rate*100);
logf(logFid, '═══════════════════\n\n');

% Per-minute timecourse
[tc_table, tc_plot] = build_group_timecourse(tc_all);
writetable(tc_table, fullfile(dataDir,'group_timecourse.csv'));
logf(logFid, 'QC/ENR: wrote group_timecourse.csv\n');

set_figure_defaults();

% Overview figs (QC) - PNG + PDF
try
    f = figure('Visible','off','Color','w','Position',[100 100 980 560]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    scatter(zeros(height(perfile_qc),1), perfile_qc.MedianDisp_postSettle_um, 36, perfile_qc.QC_PASS,'filled');
    hold on; yline(good_post_med_um,'k:');
    xlabel('(no pre-settle for GFP)'); ylabel(['median step (' umtxt ')']);
    title('Median step (whole movie; color = PASS)'); grid on;

    nexttile;
    histogram(perfile_qc.AutoSettle_min, 'BinWidth',5);
    xlabel('auto-settle time (min)'); ylabel('count'); title('Auto-settle detection (informational)'); grid on;

    nexttile;
    histogram(perfile_qc.SpikeFrac_postSettle_over_0p5um_per_s, 'BinWidth',0.005);
    xline(max_spike_frac,'k:'); xlabel(['spike fraction (>0.5 ' umtxt '/s, >=2 s runs)']);
    ylabel('count'); title('Spikes (whole movie; informational)'); grid on;

    nexttile;
    histogram(perfile_qc.Pegged_frac, 'BinWidth',0.005);
    xline(max_pegged_frac,'r:'); xlabel('pegged fraction'); ylabel('count');
    title('Near-ceiling frames (HARD gate)'); grid on;

    exportgraphics(f, fullfile(dataDir,'qc_overview.png'), 'Resolution', 300);
    exportgraphics(f, fullfile(dataDir,'qc_overview.pdf'));
    close(f);
    logf(logFid, 'Wrote qc_overview.png/pdf\n');
catch
end

% PASS/FAIL status strip - PNG + PDF
try
    [~, order] = sort(perfile_qc.QC_PASS, 'descend');
    tbl = perfile_qc(order,:);
    f2 = figure('Visible','off','Color','w','Position',[100 100 1100 540]);
    y = 1:height(tbl);
    c = double(tbl.QC_PASS);
    scatter(double(tbl.QC_PASS), y, 60, c, 'filled'); hold on;
    colormap(f2, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    set(gca,'YDir','reverse','YTick',y,'YTickLabel',tbl.File,'XTick',[0 1],'XTickLabel',{'FAIL','PASS'});
    xlim([-0.5 1.5]); ylim([0 height(tbl)+1]); grid on; box on;
    title('QC status (GFP: PASS = tracking integrity only)');
    exportgraphics(f2, fullfile(dataDir,'qc_status.png'), 'Resolution', 300);
    exportgraphics(f2, fullfile(dataDir,'qc_status.pdf'));
    close(f2);
    logf(logFid, 'Wrote qc_status.png/pdf\n');
catch
end

% Movement spectrum overview - PNG + PDF
try
    if ~isempty(psd_all)
        f_common = psd_all{1}.f(:);
        P_stack = nan(numel(f_common), numel(psd_all));
        for k=1:numel(psd_all)
            fk = psd_all{k}.f(:); Pk = psd_all{k}.P(:);
            if numel(fk)==numel(f_common) && max(abs(fk-f_common))<1e-12
                P_stack(:,k) = Pk;
            else
                P_stack(:,k) = interp1(fk,Pk,f_common,'linear','extrap');
            end
        end
        P_mean = mean(P_stack,2,'omitnan');

        f4 = figure('Visible','off','Color','w','Position',[120 120 1000 540]);
        plot(f_common, P_mean, 'LineWidth', 1.2); grid on;
        xlim([0 0.5]); xlabel('Frequency (Hz)'); ylabel('Power');
        title('Group step power spectrum (whole movie)');
        xline([band_low band_mid band_high], ':');
        exportgraphics(f4, fullfile(dataDir,'movement_spectrum_overview.png'), 'Resolution', 300);
        exportgraphics(f4, fullfile(dataDir,'movement_spectrum_overview.pdf'));
        close(f4);
        logf(logFid, 'Wrote movement_spectrum_overview.png/pdf\n');
    end
catch
end

% Movement overview montage - PNG + PDF
try
    f5 = figure('Visible','off','Color','w','Position',[100 100 1100 800]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    try
        boxchart(ones(height(perfile_qc),1), perfile_qc.MedianDisp_postSettle_um);
    catch
        boxplot(perfile_qc.MedianDisp_postSettle_um);
    end
    hold on; yline(good_post_med_um,'k:');
    ylabel(['median step (' umtxt ')']);
    title('Across movies (whole recording)');
    grid on;

    nexttile;
    try
        th = perfile_enr.DriftAngle_deg .* (pi/180);
        polarhistogram(th, 12, 'Normalization','probability'); title('Net drift direction');
    catch
        histogram(perfile_enr.DriftAngle_deg, -180:30:180); title('Net drift direction (deg)');
        xlabel('deg'); ylabel('probability');
    end

    nexttile; hold on;
    x = tc_plot.minute; y = tc_plot.mean; e = tc_plot.sem;
    if ~isempty(x)
        fill([x; flipud(x)], [y-e; flipud(y+e)], [0.9 0.9 0.9], 'EdgeColor','none');
        plot(x, y, 'k', 'LineWidth', 1.4);
        xlabel('minute'); ylabel(['per-minute median step (' umtxt ')']); title('Group timecourse'); grid on;
    end

    nexttile;
    c = double(perfile_qc.QC_PASS);
    xvals = max(perfile_qc.JitterRatio_post_cumPath_over_netDrift, 1e-6);
    scatter(xvals, perfile_qc.DriftRate_postSettle_um_per_s, 40, c,'filled');
    set(gca, 'XScale', 'log');
    xlim([1e-6, 1e3]);
    xlabel('JitterRatio (cumPath / netDrift) [log scale]');
    ylabel(['Drift rate (' umtxt '/s)']);
    title('Jitter vs drift (color = PASS)'); grid on;

    exportgraphics(f5, fullfile(dataDir,'movement_overview.png'), 'Resolution', 300);
    exportgraphics(f5, fullfile(dataDir,'movement_overview.pdf'));
    close(f5);
    logf(logFid, 'Wrote movement_overview.png/pdf\n');
catch
end

% -------- Write thresholds.tsv --------
try
    thresholds = {
        'max_pegged_frac', max_pegged_frac, 'fraction (HARD gate)';
        'peg_run_min_s', peg_run_min_s, 'seconds (HARD gate)';
        'nan_frac_max', nan_frac_max, 'fraction (HARD gate)';
        'flatline_fail', true, 'boolean (HARD gate)';
        'unstable_fail', true, 'boolean (HARD gate)';
        'good_post_med_um', good_post_med_um, 'um (soft/descriptive)';
        'spike_thr_rate_um_s', spike_thr_rate_um_s, 'um/s (soft/descriptive)';
        'max_spike_frac', max_spike_frac, 'fraction (soft/descriptive)';
        'fps_default', fps_default, 'Hz';
        'pixel_size_um', px_um_default, 'um/pixel';
    };
    thrTbl = cell2table(thresholds, 'VariableNames', {'Parameter', 'Value', 'Units_or_Description'});
    writetable(thrTbl, fullfile(dataDir, 'thresholds.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    logf(logFid, 'Wrote thresholds.tsv\n');
catch
end

% -------- Write column dictionary --------
try
    write_column_dictionary(dataDir);
    logf(logFid, 'Wrote column_dictionary_motion.csv\n');
catch
end

% -------- Write compact params JSON (provenance) --------
try
    params = struct();
    params.analyzer       = "movement_1fps_gfp_v2";
    params.created_utc    = string(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
    if isfield(versions_in,'MATLAB_release'),   params.matlab_release    = string(versions_in.MATLAB_release);   else, params.matlab_release = "n/a"; end
    if isfield(versions_in,'NoRMCorre_commit'), params.normcorre_commit  = string(versions_in.NoRMCorre_commit); else, params.normcorre_commit = "n/a"; end
    params.fps            = fps_default;
    params.px_um          = px_um_default;

    rng_state = rng;
    params.rng_type = rng_state.Type;
    params.rng_seed = rng_state.Seed;

    params.QC_thresholds_hard = struct( ...
        'max_pegged_frac',  max_pegged_frac, ...
        'peg_run_min_s',    peg_run_min_s, ...
        'nan_frac_max',     nan_frac_max, ...
        'note', 'GFP mode: only tracking integrity gates PASS');

    params.QC_thresholds_soft = struct( ...
        'good_post_med_um',     good_post_med_um, ...
        'spike_thr_rate_um_s',  spike_thr_rate_um_s, ...
        'max_spike_frac',       max_spike_frac, ...
        'imp_jump_rate_um_s',   imp_jump_rate_um_s, ...
        'imp_jump_frac_max',    imp_jump_frac_max, ...
        'note', 'Descriptive only for GFP');

    params.Settle_policy = struct( ...
        'base_settle_min', base_settle_min, ...
        'use_auto_settle_for_QC', use_auto_settle_for_QC, ...
        'note', 'GFP: analyze entire movie');

    params.PSD_bands = struct( ...
        'low',  band_low, ...
        'mid',  band_mid, ...
        'high', band_high, ...
        'note', 'Adaptive based on Nyquist frequency');

    if ~isempty(defs_in) && isstruct(defs_in) && isfield(defs_in, 'bin_width_default')
        ncdefs = struct( ...
            'bin_width',   defs_in.bin_width_default, ...
            'init_batch',  defs_in.init_batch_default, ...
            'us_fac',      defs_in.us_fac_default, ...
            'chunk_frames',defs_in.chunk_frames_default);
    else
        tmp = movement_defaults(cfg);
        ncdefs = struct('bin_width',tmp.bin_width_default,'init_batch',tmp.init_batch_default,'us_fac',tmp.us_fac_default,'chunk_frames',tmp.chunk_frames_default);
    end
    params.NoRMCorre_defaults = ncdefs;

    params.max_shift_policy = struct( ...
        'mode', cfg.max_shift_mode, ...
        'note', 'auto = 30% of FOV; rescue grows by 1.5x up to 50% FOV minus 1px', ...
        'rescue_cap_percent', 0.50, ...
        'default_window_percent', 0.30);

    params.notes = "GFP 1 Hz tracking; X/Y swap fixed Dec 2025.";

    % Record full config for reproducibility
    params.config = cfg;

    jsonStr = jsonencode(params);
    fid = fopen(fullfile(dataDir,'params_motion_1fps.json'),'w');
    if fid>0
        fwrite(fid, jsonStr, 'char'); fclose(fid);
    end
catch ME
    logf(logFid, '[WARN] Could not write params_motion_1fps.json: %s\n', ME.message);
end

end

% =====================================================================
% =========================== HELPERS =================================
% =====================================================================

function dataDirs = prompt_for_folders()
% Prompt the user to select one or more data folders interactively,
% or error with usage instructions if no GUI is available.
if usejava('desktop')
    dataDirs = {};
    while true
        folder = uigetdir('', sprintf('Select data folder %d (Cancel when done)', numel(dataDirs)+1));
        if folder == 0
            break
        end
        dataDirs{end+1} = folder; %#ok<AGROW>
    end
    if isempty(dataDirs)
        error('No folders selected. Exiting.');
    end
else
    error(['No data folders specified and no GUI available.\n' ...
           'Usage:\n' ...
           '  analyze_motion_and_QC(''/path/to/data'')\n' ...
           '  analyze_motion_and_QC({''/path/A'', ''/path/B''})\n\n' ...
           'See also: motion_analysis_config']);
end
end

function cfg = merge_with_defaults(user_cfg)
% Merge a user-provided config struct with defaults. Unknown fields
% trigger a warning. Missing fields are filled from defaults.
cfg = motion_analysis_config();
default_fields = fieldnames(cfg);
user_fields = fieldnames(user_cfg);

for i = 1:numel(user_fields)
    if ismember(user_fields{i}, default_fields)
        cfg.(user_fields{i}) = user_cfg.(user_fields{i});
    else
        warning('motion:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See motion_analysis_config() for valid fields.', ...
            user_fields{i});
    end
end
end

function logf(fid, varargin)
% Log to both console and file
fprintf(varargin{:});
if fid > 0
    fprintf(fid, varargin{:});
end
end

function set_figure_defaults()
set(groot, 'defaultAxesFontSize', 10, ...
           'defaultTextInterpreter', 'none', ...
           'defaultAxesLabelFontSizeMultiplier', 1.0, ...
           'defaultLineLineWidth', 1.1);
end

function write_column_dictionary(dataDir)
dict = {
    'motion_qc.csv', 'QC_PASS', 'boolean', 'Primary pass/fail (GFP: tracking integrity only)';
    'motion_qc.csv', 'QC_PASS_tracking', 'boolean', 'Same as QC_PASS for GFP';
    'motion_qc.csv', 'QC_soft_still', 'boolean', 'Soft check: low motion AND low spikes (informational)';
    'motion_qc.csv', 'QC_reason', '', 'Failure reason(s) or OK';
    'motion_qc.csv', 'HealthFlag', '', 'Quick status: OK, PEGGED, CEILING_RUN, etc.';
    'motion_qc.csv', 'nan_frac', 'fraction', 'Fraction of non-finite frames';
    'motion_qc.csv', 'rescue_count', 'count', 'Number of registration window rescues';
    'motion_qc.csv', 'unstable_tracking', 'boolean', 'True if rescue cap hit repeatedly';
    'motion_basic.csv', 'File', '', 'Filename of source video';
    'motion_basic.csv', 'n_frames', 'frames', 'Total frames processed';
    'motion_basic.csv', 'MedianDisp_um', 'um', 'Median per-frame step (whole movie for GFP)';
    'motion_basic.csv', 'Pegged_frac', 'fraction', 'Fraction at registration ceiling';
    'motion_basic.csv', 'Flatline_flag', 'boolean', 'True if no movement detected';
};

dictTable = cell2table(dict, 'VariableNames', {'FileName', 'Column', 'Units', 'Description'});
writetable(dictTable, fullfile(dataDir, 'column_dictionary_motion.csv'));
end

function p = Pctl(x, q)
x = x(isfinite(x));
if isempty(x), p = NaN; return; end
if exist('prctile','file')==2
    p = prctile(x, q);
else
    p = quantile(x, q/100);
end
end

function defs = movement_defaults(cfg)
% Build internal defaults struct from user config.
defs = struct();
defs.px_um               = cfg.pixel_size_um;
defs.fps                 = cfg.fps;
defs.bin_width_default   = cfg.bin_width;
defs.us_fac_default      = cfg.us_fac;
defs.init_batch_default  = cfg.init_batch;
defs.chunk_frames_default= cfg.chunk_frames;
defs.settle_min          = 10;  % not used for GFP; retained for compatibility
end

function versions = detect_versions()
try
    v = ver('MATLAB');
    if ~isempty(v) && isfield(v,'Release') && ~isempty(v.Release)
        matlabRelease = regexprep(v.Release,'[()]','');
    else
        matlabRelease = 'unknown';
    end
catch
    matlabRelease = 'unknown';
end

ncCommit = 'n/a';
try
    ncFile = which('normcorre');
    if ~isempty(ncFile)
        ncDir = fileparts(ncFile);
        [status, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', ncDir));
        if status==0
            ncCommit = strtrim(out);
        end
    end
catch
end

versions = struct('MATLAB_release', matlabRelease, 'NoRMCorre_commit', ncCommit);
end

function umtxt = get_um_label(useUnicode)
if useUnicode, umtxt = [char(181) 'm']; else, umtxt = 'um'; end
end

function max_shift_px = compute_max_shift_px(H, W, px_um, mode, val)
minHW = min(H,W);
half_lim = floor(minHW/2) - 1;
switch lower(string(mode))
    case "auto"
        max_shift_px = max(1, min(half_lim, round(0.30 * minHW)));
    case "auto_conservative"
        max_shift_px = max(1, min(half_lim, round(50 / px_um)));
    case "auto_loose"
        max_shift_px = half_lim;
    case "fixed_um"
        if isempty(val), error('max_shift_value (um) required for fixed_um mode.'); end
        v = val(:)'; if numel(v)==1, v=[v v]; end
        ms = round(v./px_um);
        max_shift_px = [max(1, min(ms(1), half_lim)), max(1, min(ms(2), half_lim))];
        if numel(unique(max_shift_px)) == 1
            max_shift_px = max_shift_px(1);
        end
    case "fixed_px"
        if isempty(val), error('max_shift_value (px) required for fixed_px mode.'); end
        v = val(:)'; if numel(v)==1, v=[v v]; end
        ms = round(v);
        max_shift_px = [max(1, min(ms(1), half_lim)), max(1, min(ms(2), half_lim))];
        if numel(unique(max_shift_px)) == 1
            max_shift_px = max_shift_px(1);
        end
    otherwise
        error('Unknown max_shift_mode: %s', mode);
end
end

function n_bursts = burst_count(x_um, thr_um, fps)
if nargin<3 || isempty(x_um) || isempty(thr_um), n_bursts = NaN; return; end
x = x_um > thr_um; rises = [false; x(2:end) & ~x(1:end-1)];
n_bursts = sum(rises);
end

function rate_per_min = burst_rate(x_um, thr_um, fps)
if nargin<3 || isempty(x_um) || isempty(thr_um) || isempty(fps)
    rate_per_min = NaN; return;
end
x = x_um > thr_um; rises = [false; x(2:end) & ~x(1:end-1)];
n_bursts = sum(rises); dur_min = max(1e-9, numel(x_um)/ (60*fps));
rate_per_min = n_bursts / dur_min;
end

function [f, Pxx] = compute_psd_safe(x, fps)
x = x(:); x = x - mean(x,'omitnan'); x(~isfinite(x)) = 0;
if exist('pwelch','file')==2
    try
        [Pxx_out, f_out] = pwelch(x, hamming(min(128, max(16, floor(numel(x)/4)))), [], [], fps);
        f = f_out; Pxx = Pxx_out;
        return;
    catch
    end
end
N = 2^nextpow2(max(64, numel(x)));
X = fft(x, N); P = (abs(X).^2)/N; f = (0:N-1)'*(fps/N);
half = floor(N/2)+1; f = f(1:half); Pxx = P(1:half);
end

function [f0, p0] = dominant_peak(f, Pxx, frange)
if isempty(f), f0=0; p0=0; return; end
mask = (f>=frange(1) & f<=frange(2));
if ~any(mask), f0=0; p0=0; return; end
[~,idx] = max(Pxx(mask)); ff = f(mask); pp = Pxx(mask);
f0 = ff(idx); p0 = pp(idx);
end

function frac = band_power_frac(f, Pxx, band)
if isempty(f) || all(~isfinite(Pxx)), frac = NaN; return; end
Ptot = trapz(f, Pxx); if Ptot<=0, frac = NaN; return; end
m = (f>=band(1) & f<=band(2)); Pb = trapz(f(m), Pxx(m));
frac = Pb / Ptot;
end

function frac = frac_above_threshold(x, thr, min_consec)
x = x(:) > thr;
if ~isfinite(min_consec) || min_consec <= 1
    frac = mean(x); return;
end
d = diff([0; x; 0]);
starts = find(d==1); ends = find(d==-1) - 1;
keep = false(size(x));
for k = 1:numel(starts)
    if (ends(k)-starts(k)+1) >= min_consec
        keep(starts(k):ends(k)) = true;
    end
end
frac = mean(keep);
end

function [minute_idx, per_min_median] = per_minute_median(steps_um, fps)
nF = numel(steps_um);
mins = floor(((0:nF-1).')/(60*fps)) + 1;
u = unique(mins);
per_min_median = zeros(numel(u),1);
for ii=1:numel(u)
    per_min_median(ii) = median(steps_um(mins==u(ii)), 'omitnan');
end
minute_idx = u;
end

function [tc_table, tc_plot] = build_group_timecourse(tc_all)
maxM = 0; for ii=1:numel(tc_all), maxM = max(maxM, max(tc_all{ii}.minute)); end
M = (1:maxM).';
vals = NaN(maxM, numel(tc_all));
for j=1:numel(tc_all), m = tc_all{j}.minute; v = tc_all{j}.file_median; vals(m,j) = v; end
mean_v = mean(vals, 2, 'omitnan'); n_v = sum(~isnan(vals),2);
sem_v  = std(vals, 0, 2, 'omitnan') ./ sqrt(max(1,n_v));
tc_table = table(M, n_v, mean_v, sem_v, 'VariableNames', ...
    {'minute','n_files_active','mean_perMinuteMedian_um','sem_perMinuteMedian_um'});
tc_plot = struct('minute',M,'mean',mean_v,'sem',sem_v);
end

function [dx, dy] = extract_dx_dy(shifts)
% Extract dx (horizontal/X) and dy (vertical/Y) from NoRMCorre shifts
% IMPORTANT: NoRMCorre uses [row, col] = [Y, X] convention
%
% OUTPUT CONVENTION (MOTION direction):
%   dx > 0 = sample moved RIGHT
%   dy > 0 = sample moved DOWN
%
% NoRMCorre has asymmetric sign conventions:
%   - X output is inverted, so we NEGATE
%   - Y output is already motion direction, keep as-is
dx = []; dy = [];
if isstruct(shifts)
    if isfield(shifts,'x') && isfield(shifts,'y')
        dx = -double(shifts.x(:));
        dy = double(shifts.y(:));
        return;
    end
    if isfield(shifts,'shifts')
        try
            tmp = arrayfun(@(s) s.shifts(:), shifts, 'UniformOutput', false);
            tmp = cat(2, tmp{:});
            dx = -double(tmp(2,:)).';
            dy = double(tmp(1,:)).';
            return;
        catch
        end
    end
    fns = fieldnames(shifts);
    for ii = 1:numel(fns)
        [dx,dy] = coerce_numeric_to_dxdy(shifts.(fns{ii}));
        if ~isempty(dx), return; end
    end
elseif isnumeric(shifts)
    [dx,dy] = coerce_numeric_to_dxdy(shifts);
end
end

function [dx, dy] = coerce_numeric_to_dxdy(val)
% Coerce numeric array to dx, dy
% NoRMCorre convention: [row, col] = [Y, X]
% X is negated (inverted in NoRMCorre), Y kept as-is
dx = []; dy = [];
if ~isnumeric(val), return; end
sz = size(val);
if numel(sz)==2 && sz(2)==2
    dx = -double(val(:,2));
    dy = double(val(:,1));
elseif numel(sz)==2 && sz(1)==2
    dx = -double(val(2,:)).';
    dy = double(val(1,:)).';
elseif isvector(val)
    dx = -double(val(:));     dy = zeros(size(dx));
end
end

function fr = readFrame_gray(vr)
fr = readFrame(vr);
if ndims(fr)==3
    if exist('rgb2gray','file')==2, fr = rgb2gray(fr); else, fr = fr(:,:,2); end
end
end

function A = rescale01(A)
A = single(A); mn = min(A(:)); mx = max(A(:)); rng_val = mx - mn;
if rng_val <= 0, A = zeros(size(A), 'like', A); else, A = (A - mn) ./ rng_val; end
end

function S = compute_summary_from_csv(csvPath, pxSize_um, fps_default, max_shift_px_in)
T = readtable(csvPath);
if ismember('px_um', T.Properties.VariableNames),       pxSize_um   = T.px_um(1);     end
if ismember('fps_used', T.Properties.VariableNames),    fps_default = T.fps_used(1);  end

if ismember('max_shift_um_initial', T.Properties.VariableNames)
    maxShift_um_initial = T.max_shift_um_initial(1);
elseif ismember('max_shift_px_initial', T.Properties.VariableNames)
    maxShift_um_initial = T.max_shift_px_initial(1) * pxSize_um;
elseif ~isnan(max_shift_px_in)
    maxShift_um_initial = max_shift_px_in * pxSize_um;
else
    maxShift_um_initial = NaN;
end

if ismember('max_shift_um', T.Properties.VariableNames)
    maxShift_um_effective = T.max_shift_um(1);
elseif ismember('max_shift_px', T.Properties.VariableNames)
    maxShift_um_effective = T.max_shift_px(1) * pxSize_um;
else
    maxShift_um_effective = maxShift_um_initial;
end

if all(ismember({'dx_um','dy_um'}, T.Properties.VariableNames))
    dx_um = T.dx_um; dy_um = T.dy_um;
else
    dx_um = T.dx_px * pxSize_um; dy_um = T.dy_px * pxSize_um;
end
steps_um = hypot([0; diff(dx_um)], [0; diff(dy_um)]);
nF = height(T);

win  = max(3, 2*floor((60*fps_default)/2)+1);
try
    dx_s = medfilt1(dx_um, win, 'omitnan', 'truncate');
    dy_s = medfilt1(dy_um, win, 'omitnan', 'truncate');
catch
    dx_s = movmedian(dx_um, win, 'omitnan');
    dy_s = movmedian(dy_um, win, 'omitnan');
end
net_drift  = hypot(dx_s(end)-dx_s(1), dy_s(end)-dy_s(1));
duration_s = (nF - 1) / max(1, fps_default);
drift_rate = net_drift / max(1e-9, duration_s);

med_step  = median(steps_um, 'omitnan');
p95_step  = Pctl(steps_um, 95);
cum_path  = sum(steps_um, 'omitnan');

settle_min_gfp = 0;
med_step_ps   = med_step;
p95_step_ps   = p95_step;
cum_path_ps   = cum_path;
net_drift_ps  = net_drift;
drift_rate_ps = drift_rate;

[~, base] = fileparts(csvPath);
fileStem = regexprep(base, '_shifts$', '');

S = { [fileStem '.avi'], nF, ...  % NOTE: extension is cosmetic for summary table label
    med_step, p95_step, cum_path, net_drift, drift_rate, NaN, NaN, ...
    med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, NaN, NaN, ...
    fps_default, pxSize_um, ...
    maxShift_um_initial, maxShift_um_effective, ...
    settle_min_gfp };
end

function v = steps_post_if(post_idx, steps_um)
if isempty(post_idx) || numel(post_idx)<2, v = []; else, v = steps_um(post_idx); end
end

function m = max_run(xlog)
x = logical(xlog(:));
d = diff([false; x; false]);
starts = find(d==1); ends = find(d==-1)-1;
if isempty(starts), m = 0; else, m = max(ends - starts + 1); end
end

function [starts, ends] = find_runs(x)
x = x(:) ~= 0;
d = diff([false; x; false]);
starts = find(d==1);
ends = find(d==-1)-1;
end

function shade_segments(t, mask, rgb)
[starts, ends] = find_runs(mask);
yl = ylim;
for ii=1:numel(starts)
    if starts(ii) <= numel(t) && ends(ii) <= numel(t)
        xs = [t(starts(ii)) t(ends(ii)) t(ends(ii)) t(starts(ii))];
        ys = [yl(1) yl(1) yl(2) yl(2)];
        patch(xs, ys, rgb, 'EdgeColor','none', 'FaceAlpha',0.25);
    end
end
lines = findobj(gca,'Type','line');
if ~isempty(lines)
    uistack(lines,'top');
end
end