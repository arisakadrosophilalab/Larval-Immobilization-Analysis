function analyze_motion_and_QC_neuro(dataDirs, varargin)
% ANALYZE_MOTION_AND_QC_NEURO  Rigid motion tracking and dual quality
%   control for neurological calcium imaging recordings (~33.33 Hz).
%
%   Uses NoRMCorre rigid registration to estimate frame-by-frame
%   translational shifts, computes motion metrics, and applies a dual
%   QC system (MOVEMENT_VALID + ROI_ELIGIBLE).
%
%   USAGE:
%     analyze_motion_and_QC_neuro()                    % Opens folder picker
%     analyze_motion_and_QC_neuro('/path/to/videos')   % Process one folder
%     analyze_motion_and_QC_neuro({'/A', '/B'})         % Process multiple folders
%
%     cfg = neuro_motion_config();                     % Get default config
%     cfg.pixel_size_um = 0.5;                         % Modify as needed
%     cfg.fps = 30;
%     analyze_motion_and_QC_neuro('/path', cfg)         % Use custom config
%
%   INPUTS:
%     dataDirs - (optional) Char, string, or cell array of folder paths
%                containing video files. If omitted, a folder picker opens.
%     cfg      - (optional) Configuration struct from neuro_motion_config().
%                If omitted, defaults are used. Unknown fields trigger a warning.
%
%   OUTPUTS (written to each data folder):
%     Per-file outputs:
%       *_shifts.csv             Per-frame shift estimates
%       *_motion_trace.png/pdf   Motion trace plots
%       *_rescue_meta.json       Registration rescue event log
%       *_params.json            Effective per-file registration parameters
%       *_pegged_runs.csv        Contiguous ceiling-hit episodes
%       *_reposition_mask.csv    Repositioning event mask (if detected)
%
%     Per-folder outputs:
%       motion_summary.csv       Summary statistics across files
%       motion_basic.csv         Basic motion metrics
%       motion_qc.csv            Dual QC flags and metrics
%       motion_enriched.csv      Extended metrics (PSD, bursts, etc.)
%       motion_group_summary.csv Group-level aggregates
%       group_timecourse.csv     Per-minute motion timecourse
%       roi_candidates.csv/json  Pipeline handoff with dual QC flags
%       movement_valid_list.txt  Files passing MOVEMENT_VALID
%       roi_eligible_list.txt    Files passing ROI_ELIGIBLE
%       params_motion_33fps.json Full parameter provenance record
%       thresholds.tsv           QC threshold documentation
%       column_dictionary_motion.csv  Column definitions
%       qc_dual_status.png/pdf   Dual QC status strip
%       qc_boundary_excursion.png/pdf  Boundary excursion plots
%       qc_overview.png/pdf      QC summary plots
%       movement_overview.png/pdf  Motion overview montage
%       movement_spectrum_overview.png/pdf  Group PSD
%       motion_run.log           Processing log
%
%   DUAL QC SYSTEM:
%     MOVEMENT_VALID - Can we trust motion measurements for comparing
%       immobilization techniques? Checks: nan_frac, pegged_frac,
%       flatline, unstable_tracking. Does NOT check boundary excursion.
%     ROI_ELIGIBLE - Is this file suitable for ROI detection and calcium
%       trace extraction? Requires MOVEMENT_VALID plus brain stayed in FOV.
%     A file can have MOVEMENT_VALID=true but ROI_ELIGIBLE=false if the
%     brain partially exited the FOV.
%
%   CONFIGURATION:
%     All parameters (pixel size, frame rate, file pattern, QC thresholds,
%     NoRMCorre settings, boundary checks, settle policy) are set via the
%     config struct. See neuro_motion_config() for the full list.
%
%   DEPENDENCIES:
%     - MATLAB R2020b or later (recommended: R2024b)
%     - NoRMCorre (https://github.com/flatironinstitute/NoRMCorre)
%     - Statistics and Machine Learning Toolbox (optional; for prctile)
%     - Signal Processing Toolbox (optional; for medfilt1, pwelch)
%
%   SHIFT CONVENTION (CSV output):
%     dx_px > 0 = sample moved RIGHT
%     dy_px > 0 = sample moved UP
%     This is MOTION direction, intuitive for visualization.
%     The pipeline applies [-dx, +dy] for correction (Cartesian to image coords).
%
%   NOTE: AVI containers typically truncate fps to integer (33 Hz instead
%   of 33.33 Hz). This pipeline uses fps_override=33.33 by default.
%   10,000 frames at 33.33 Hz = 300 seconds = 5.00 minutes exactly.
%
%   See also: neuro_motion_config

% ---------- Parse config ----------
if nargin >= 2 && isstruct(varargin{1})
    cfg = merge_with_defaults(varargin{1});
elseif nargin >= 2
    error(['Unexpected second argument (type: %s). Pass a config struct from\n' ...
           'neuro_motion_config(), or call with just the data folder path.'], ...
           class(varargin{1}));
else
    cfg = neuro_motion_config();
end

% ---------- Fail-fast toolbox checks ----------
assert(exist('VideoReader','class') == 8, ...
    'VideoReader not found. A standard MATLAB installation is required.');
if exist('pwelch','file') ~= 2
    warning('Signal Processing Toolbox (pwelch) not found. Using FFT fallback for PSD.');
end

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

% normalize to cellstr of char paths
if isstring(dataDirs), dataDirs = cellstr(dataDirs); end
if ~iscell(dataDirs),  dataDirs  = {dataDirs};       end
dataDirs = cellfun(@(p) char(string(p)), dataDirs, 'UniformOutput', false);

% filter non-existing with warning (but continue others)
ok = true(size(dataDirs));
for i=1:numel(dataDirs)
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
        if k==6 && numel(files) > 6
            fprintf('   ... (%d more)\n', numel(files)-5);
            break
        else
            fprintf('   %s\n', fullfile(files(k).folder, files(k).name));
        end
    end
    fprintf('\n=== Processing folder [%d/%d]: %s (%d files) ===\n\n', ...
        i, N, dataDir, numel(files));

    % analyze + QC for this folder
    analyze_motion_and_QC_neuro_core(dataDir, cfg);
end

fprintf('\nAll done: analyzer + QC + enriched complete for %d folder(s).\n', N);
end

% =====================================================================
% =============== CORE ENGINE (with logging + improvements) ===========
% =====================================================================
function analyze_motion_and_QC_neuro_core(dataDir, cfg)
% One-shot: run rigid NoRMCorre + QC + enriched on video files in dataDir.

% -------- Open log file --------
logPath = fullfile(dataDir, 'motion_run.log');
logFid = fopen(logPath, 'a');
if logFid < 0
    warning('Could not open log file: %s', logPath);
    logFid = -1;
end

logf(logFid, '\n════════════════════════════════════════════════════════════════\n');
logf(logFid, 'Neuro Motion Analysis Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
logf(logFid, 'Folder: %s\n', dataDir);
logf(logFid, 'Config: pixel_size_um=%.4f, fps=%g, fps_override=%s, file_pattern=%s\n', ...
    cfg.pixel_size_um, cfg.fps, mat2str(cfg.fps_override), cfg.file_pattern);
logf(logFid, '════════════════════════════════════════════════════════════════\n\n');

try
    % ---------- run analyzer ----------
    [defaults_struct, versions_struct] = run_analyzer_neuro(dataDir, cfg, logFid);

    % ---------- QC + ENRICHED ----------
    qc_and_enriched_motion_folder_neuro(dataDir, cfg, defaults_struct, versions_struct, logFid);
    
    logf(logFid, '\n════════════════════════════════════════════════════════════════\n');
    logf(logFid, 'Neuro Motion Analysis Completed: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
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
% ===================== ANALYZER (NEURO) ==============================
% =====================================================================
function [defaults_out, versions_out] = run_analyzer_neuro(dataDir, cfg, logFid)
% Analyzer function: NoRMCorre registration for neuro recordings.

% Build defaults from config
defs = neuro_defaults(cfg);
pxSize_um = defs.px_um;
fps       = defs.fps;

% Report defaults outward for JSON provenance
defaults_out = defs;
defaults_out.fps_override = cfg.fps_override;

% Unpack config fields used in this function
fps_override     = cfg.fps_override;
spike_thr_rate_um_s = cfg.qc_spike_thr_um_s;
spike_min_consec_fr = cfg.qc_spike_min_consec_fr;

% for labels
umtxt = get_um_label(cfg.use_unicode_labels);

% Check NoRMCorre is available
assert(exist('NoRMCorreSetParms','file')==2 && exist('normcorre','file')==2, ...
    'NoRMCorre not found on path. Add it first.');

% versions (MATLAB + NoRMCorre git short)
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
    
    % Guard VideoReader creation for corrupt/weird AVI
    try
        vr = VideoReader(aviPath);
    catch ME
        logf(logFid, '[WARN] VideoReader failed for %s: %s\n', aviPath, ME.message);
        continue
    end
    
    % Determine effective fps: use override if provided, else measured
    fps_meas = vr.FrameRate;
    if ~isempty(fps_override) && fps_override > 0
        fps_eff = fps_override;
        logf(logFid, '  Using fps_override=%.4f (AVI reports %.1f)\n', fps_eff, fps_meas);
    else
        fps_eff = fps_meas;
    end
    fps_exp = fps;

    % reuse if present
    if cfg.skip_if_present && exist(outCSV,'file')==2
        logf(logFid, 'Found existing shifts CSV; skipping registration: %s\n', outCSV);
        S = compute_summary_from_csv(outCSV, pxSize_um, fps, NaN, defs.settle_min);
        if numel(S) ~= 20
            error('Summary row for %s has %d columns (expected 20).', files(f).name, numel(S));
        end
        allSummary(end+1,:) = S; %#ok<AGROW>
        continue;
    end

    H = vr.Height; W = vr.Width;
    nF_est = floor(vr.Duration * max(1, vr.FrameRate));
    logf(logFid, 'Estimated frames: %d; size %dx%d\n', nF_est, H, W);
    
    % Short-clip awareness at 33 Hz
    is_short_clip = nF_est < round(10 * fps_eff);
    if is_short_clip
        logf(logFid, '  Short clip detected (<10 seconds), adjusting parameters...\n');
        this_bin_width = max(200, round(defs.bin_width_default/2));
        this_init_batch = max(200, round(defs.init_batch_default/2));
    else
        this_bin_width = defs.bin_width_default;
        this_init_batch = defs.init_batch_default;
    end

    % max shift with anisotropic support + rescue tracking
    max_shift_px_initial = compute_max_shift_px(H, W, pxSize_um, cfg.max_shift_mode, cfg.max_shift_value);
    rescue_cap_px = compute_rescue_cap(H, W, max_shift_px_initial);
    
    max_shift_px = max_shift_px_initial;
    if isscalar(max_shift_px)
        maxShift_um_scalar = max_shift_px * pxSize_um;
    else
        maxShift_um_scalar = max(max_shift_px) * pxSize_um;
    end
    
    fmt_ms = @(v) mat2str(v);
    logf(logFid, '  Using max_shift=%s px (~%s %s) [rescue cap=%s px]\n', ...
        fmt_ms(max_shift_px), fmt_ms(max_shift_px.*pxSize_um), umtxt, fmt_ms(rescue_cap_px));
    
    rescue_count = 0;
    rescue_history = {};
    unstable_tracking = false;

    % RAM caps
    bytesPerFrame   = H * W * 4;
    maxChunkByRAM   = max(120, floor((cfg.target_gb*1e9) / bytesPerFrame));
    chunk_frames    = min(defs.chunk_frames_default, maxChunkByRAM);
    this_init_batch = min(this_init_batch, chunk_frames);
    if this_init_batch < defs.init_batch_default
        logf(logFid, '  Capping init_batch from %d to %d for RAM budget (~%.1f GB target)\n', defs.init_batch_default, this_init_batch, cfg.target_gb);
    end
    if chunk_frames < defs.chunk_frames_default
        logf(logFid, '  Capping chunk_frames from %d to %d for RAM budget (~%.1f GB target)\n', defs.chunk_frames_default, chunk_frames, cfg.target_gb);
    end

    % NoRMCorre opts
    options_mc = NoRMCorreSetParms('d1',H,'d2',W, ...
        'bin_width',this_bin_width, ...
        'max_shift',max_shift_px, ...
        'us_fac',defs.us_fac_default, ...
        'init_batch',this_init_batch);
    
    us_fac_template_used = options_mc.us_fac;
    us_fac_stream_used = options_mc.us_fac;

    % Build template from quietest ~10s window
    win_s = 10; 
    winN = max(200, round(win_s*fps_eff));
    vr.CurrentTime = 0;
    
    bytesPerFrame = H*W*4;
    maxInitByRAM = max(200, floor((cfg.target_gb*1e9)/(bytesPerFrame*3)));
    bufN = min(max(this_init_batch + round(10*fps_eff), 800), round(vr.FrameRate*30));
    bufN = min(bufN, maxInitByRAM);
    
    B = zeros(H,W,bufN,'single'); 
    n=0;
    while hasFrame(vr) && n<bufN
        n=n+1; 
        B(:,:,n)=single(readFrame_gray(vr)); 
    end
    B = B(:,:,1:n); 
    if n<winN, winN = n; end

    d = [zeros(1,1,'single'), squeeze(median(abs(diff(reshape(B,[],n),1,2)),1))];
    s = movmedian(d, round(0.5*fps_eff));
    w = winN; 
    mm = movmean(s, w);
    [~,i0] = min(mm(1:n-w+1)); 
    i1 = i0+w-1;
    initBuf = B(:,:,i0:i1); 
    clear B d s mm;
    logf(logFid, 'Building initial template from quiet window [frames %d-%d] (~%.1f s)\n', i0, i1, w/fps_eff);

    template = []; shifts_init = [];
    try
        [~, shifts_init, template] = normcorre(initBuf, options_mc); %#ok<ASGLU>
    catch ME
        if contains(ME.message,'Out of memory','IgnoreCase',true)
            logf(logFid, '  OOM during template; retry with lighter settings (us_fac/2, init_batch/2)...\n');
            this_us_fac = max(10, round(defs.us_fac_default/2));
            us_fac_template_used = this_us_fac;
            this_init_batch = max(120, round(this_init_batch/2));
            options_mc = NoRMCorreSetParms('d1',H,'d2',W, ...
                'bin_width',this_bin_width,'max_shift',max_shift_px, ...
                'us_fac',this_us_fac,'init_batch',this_init_batch);
            initBuf = initBuf(:,:,1:min(size(initBuf,3), this_init_batch));
            [~, shifts_init, template] = normcorre(initBuf, options_mc); %#ok<ASGLU>
        else
            rethrow(ME);
        end
    end
    clear initBuf;
    
    % PREALLOCATE dx/dy arrays
    dx_all = zeros(nF_est, 1);
    dy_all = zeros(nF_est, 1);
    writeIdx = 0;
    
    vr.CurrentTime = 0;

    options_stream = options_mc; options_stream.init_batch = 0;
    chunk = zeros(H, W, chunk_frames, 'single'); nIn = 0;

    while hasFrame(vr)
        fr = readFrame_gray(vr);
        fr = single(fr);
        nIn = nIn + 1; chunk(:,:,nIn) = fr;

        if nIn == chunk_frames || ~hasFrame(vr)
            Yc = chunk(:,:,1:nIn);

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
                        us_fac_stream_used = options_stream.us_fac;
                        logf(logFid, '  [OOM] retry with lower us_fac=%d\n', options_stream.us_fac);
                    else
                        rethrow(ME);
                    end
                end
            end

            % Bounded near-ceiling rescue
            [dx_try, dy_try] = extract_dx_dy(shifts_c);
            chunk_rescue_count = 0;
            
            while isfield(options_stream,'max_shift') && ~isempty(dx_try) && chunk_rescue_count < 5
                ms = options_stream.max_shift;
                if isscalar(ms), msy = ms; msx = ms; else, msy = ms(1); msx = ms(2); end
                if isscalar(rescue_cap_px), capy = rescue_cap_px; capx = rescue_cap_px; 
                else, capy = rescue_cap_px(1); capx = rescue_cap_px(2); end
                
                hit = any(abs(dx_try) > 0.90*msx) || any(abs(dy_try) > 0.90*msy);
                if ~hit, break; end
                
                bump_y = min(capy, max(msy+10, 2*msy));
                bump_x = min(capx, max(msx+10, 2*msx));
                new_ms_vec = [bump_y, bump_x];
                
                if all([bump_y <= msy, bump_x <= msx]) || any([bump_y > capy, bump_x > capx])
                    break;
                end
                
                try
                    opts_bump = options_stream;
                    if isscalar(ms)
                        opts_bump.max_shift = max(new_ms_vec);
                    else
                        opts_bump.max_shift = new_ms_vec;
                    end
                    logf(logFid, '  [RESCUE] near ceiling %s -> max_shift=%s px\n', ...
                        fmt_ms([msy msx]), fmt_ms(opts_bump.max_shift));
                    [~, shifts_c2, template] = normcorre(Yc, opts_bump, template); %#ok<ASGLU>
                    shifts_c = shifts_c2;
                    options_stream.max_shift = opts_bump.max_shift;
                    
                    rescue_count = rescue_count + 1;
                    chunk_rescue_count = chunk_rescue_count + 1;
                    frame_start = writeIdx + 1;
                    frame_end = frame_start + size(Yc,3) - 1;
                    rescue_history{end+1} = sprintf('%s->%s px, frames %d-%d', ...
                        fmt_ms([msy msx]), fmt_ms(opts_bump.max_shift), frame_start, frame_end); %#ok<AGROW>
                    [dx_try, dy_try] = extract_dx_dy(shifts_c);
                catch
                    break
                end
            end
            
            if chunk_rescue_count >= 5
                ms = options_stream.max_shift;
                if isscalar(ms), msy = ms; msx = ms; else, msy = ms(1); msx = ms(2); end
                if any(abs(dx_try) > 0.90*msx) || any(abs(dy_try) > 0.90*msy)
                    logf(logFid, '  [RESCUE] Hit ceiling after %d bumps, marking UNSTABLE\n', chunk_rescue_count);
                    unstable_tracking = true;
                end
            end

            [dx_c, dy_c] = extract_dx_dy(shifts_c);
            nNew = numel(dx_c);
            
            % Write to preallocated arrays, expand if needed
            if writeIdx + nNew > numel(dx_all)
                dx_all = [dx_all; zeros(nNew + 1000, 1)]; %#ok<AGROW>
                dy_all = [dy_all; zeros(nNew + 1000, 1)]; %#ok<AGROW>
            end
            dx_all(writeIdx+1 : writeIdx+nNew) = dx_c(:);
            dy_all(writeIdx+1 : writeIdx+nNew) = dy_c(:);
            writeIdx = writeIdx + nNew;
            
            nIn = 0;
        end
    end

    % Truncate to actual size
    dx_all = dx_all(1:writeIdx);
    dy_all = dy_all(1:writeIdx);
    nF = numel(dx_all);
    
    if abs(nF - nF_est) / max(1, nF_est) > 0.1
        logf(logFid, '[WARN] Frame count differs from duration estimate by >10%% (actual=%d, est=%d)\n', nF, nF_est);
    end
    
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
    
    logf(logFid, 'Total frames processed: %d (us_fac_template=%d, us_fac_stream=%d, init_batch=%d, chunk<=%d)\n', ...
        nF, us_fac_template_used, us_fac_stream_used, this_init_batch, chunk_frames);
    if unstable_tracking, unstable_tag = ' [UNSTABLE]'; else, unstable_tag = ''; end
    logf(logFid, '  Max shift: initial=%s px, effective=%s px, rescues=%d%s\n\n', ...
        fmt_ms(max_shift_px_initial), fmt_ms(eff_ms), rescue_count, unstable_tag);

    % save per-frame CSV
    dx_um = dx_all * pxSize_um; dy_um = dy_all * pxSize_um;
    steps_um = hypot([0; diff(dx_um)], [0; diff(dy_um)]);

    t_s   = (0:nF-1).'/fps_eff;
    t_min = t_s/60;
    
    if isscalar(max_shift_px_initial)
        init_max_shift_um = max_shift_px_initial * pxSize_um;
    else
        init_max_shift_um = max(max_shift_px_initial) * pxSize_um;
    end
    
    T = table((1:nF)', t_s, t_min, dx_all, dy_all, dx_um, dy_um, steps_um, ...
        repmat(max_shift_px_initial, nF,1), ...
        repmat(init_max_shift_um, nF,1), ...
        repmat(eff_max_shift_px, nF,1), repmat(eff_maxShift_um, nF,1), ...
        repmat(rescue_count, nF,1), ...
        repmat(logical(unstable_tracking), nF,1), ...
        repmat(fps_meas, nF,1), repmat(fps_eff, nF,1), repmat(fps_exp, nF,1), repmat(pxSize_um, nF,1), ...
        repmat(us_fac_template_used, nF,1), repmat(us_fac_stream_used, nF,1), ...
        'VariableNames', {'Frame','time_s','time_min','dx_px','dy_px','dx_um','dy_um','displacement_um', ...
        'max_shift_px_initial','max_shift_um_initial','max_shift_px','max_shift_um', ...
        'rescue_count','unstable_tracking', ...
        'fps_measured','fps_used','fps_expected','px_um', ...
        'us_fac_template_used','us_fac_stream_used'});
    
    if exist('eff_ms','var') && ~isscalar(eff_ms)
        T.max_shift_px_y = repmat(eff_ms(1), nF, 1);
        T.max_shift_px_x = repmat(eff_ms(2), nF, 1);
        T.max_shift_um_y = repmat(eff_ms(1) * pxSize_um, nF, 1);
        T.max_shift_um_x = repmat(eff_ms(2) * pxSize_um, nF, 1);
    end
    
    T.MATLAB_release   = repmat(string(versions_out.MATLAB_release), nF, 1);
    T.NoRMCorre_commit = repmat(string(versions_out.NoRMCorre_commit), nF, 1);
    
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
    
    % Write per-file parameters JSON
    file_params = struct('bin_width', this_bin_width, ...
                        'init_batch', this_init_batch, ...
                        'us_fac_default', defs.us_fac_default, ...
                        'us_fac_template_used', us_fac_template_used, ...
                        'us_fac_stream_used', us_fac_stream_used, ...
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
        peggedTable = table(peg_starts, peg_ends, (peg_ends-peg_starts+1)/fps_eff, ...
            'VariableNames', {'start_idx','end_idx','duration_s'});
        writetable(peggedTable, fullfile(files(f).folder, [stem '_pegged_runs.csv']));
    end

    % quick plot with shading (PNG + PDF)
    try
        set_figure_defaults();
        
        lim_dx   = Pctl(abs(dx_um), 99);
        lim_dy   = Pctl(abs(dy_um), 99);
        lim_step = Pctl(steps_um, 99);
        pad = 1.15;
        ydx   = max(lim_dx*pad, 5);
        ydy   = max(lim_dy*pad, 5);
        ystep = max(lim_step*pad, 0.5);
        
        spike_thr_um = spike_thr_rate_um_s / fps_eff;
        spike_rate_used = spike_thr_rate_um_s;
        if isempty(spike_min_consec_fr)
            spike_min_consec_eff = max(6, round(0.18*fps_eff));
        else
            spike_min_consec_eff = spike_min_consec_fr;
        end
        is_spike = false(size(steps_um));
        x_spike = steps_um > spike_thr_um;
        d_spike = diff([0; x_spike; 0]);
        spike_starts = find(d_spike==1); spike_ends = find(d_spike==-1)-1;
        for kk = 1:numel(spike_starts)
            if (spike_ends(kk)-spike_starts(kk)+1) >= spike_min_consec_eff
                is_spike(spike_starts(kk):spike_ends(kk)) = true;
            end
        end

        fign = figure('Visible','off','Color','w','Position',[200 200 980 760]);
        tEnd = t_min(end);
        
        subplot(3,1,1); plot(t_min, dx_um, 'LineWidth', 1); grid on;
        ylabel(['Delta x (' umtxt ')']); title(strrep(files(f).name,'_','\_'));
        ylim([-ydx ydx]); xlim([0, tEnd]);

        subplot(3,1,2); plot(t_min, dy_um, 'LineWidth', 1); grid on;
        ylabel(['Delta y (' umtxt ')']); 
        ylim([-ydy ydy]); xlim([0, tEnd]);

        subplot(3,1,3); plot(t_min, steps_um, 'LineWidth', 1); grid on; hold on;
        ylabel(['step (' umtxt ')']); xlabel('time (min)'); 
        ylim([0 ystep]); xlim([0, tEnd]);
        shade_segments(t_min, is_pegged, [1 0.85 0.85]);
        shade_segments(t_min, is_spike,  [1 0.95 0.85]);
        
        text(0.01*tEnd, 0.95*ystep, sprintf('spikes: >%.1f %s/s, ≥%d fr', ...
            spike_rate_used, umtxt, spike_min_consec_eff), ...
            'FontSize',8, 'VerticalAlignment','top');

        if defs.settle_min <= tEnd
            xline(defs.settle_min,':','Settle','LabelHorizontalAlignment','left');
        end

        outPNG = fullfile(files(f).folder, sprintf('%s_motion_trace.png', stem));
        outPDF = fullfile(files(f).folder, sprintf('%s_motion_trace.pdf', stem));
        exportgraphics(fign, outPNG, 'Resolution', 300);
        exportgraphics(fign, outPDF);
        close(fign);
    catch
    end

    % summary row computation
    win = max(3, 2*floor((60*fps_eff)/2)+1);
    try
        dx_s = medfilt1(dx_um, win, 'omitnan', 'truncate');
        dy_s = medfilt1(dy_um, win, 'omitnan', 'truncate');
    catch
        dx_s = movmedian(dx_um, win, 'omitnan');
        dy_s = movmedian(dy_um, win, 'omitnan');
    end
    net_drift  = hypot(dx_s(end)-dx_s(1), dy_s(end)-dy_s(1));
    duration_s = (nF - 1) / max(1, fps_eff);
    drift_rate = net_drift / max(1e-9, duration_s);

    med_step  = median(steps_um, 'omitnan');
    p95_step  = Pctl(steps_um,95);
    cum_path  = sum(steps_um,'omitnan');

    settle_frames = min(nF, round(defs.settle_min * 60 * fps_eff));
    idx = (settle_frames+1):nF;
    if numel(idx) >= 2
        med_step_ps   = median(steps_um(idx), 'omitnan');
        p95_step_ps   = Pctl(steps_um(idx), 95);
        cum_path_ps   = sum(steps_um(idx), 'omitnan');
        net_drift_ps  = hypot(dx_s(idx(end))-dx_s(idx(1)), dy_s(idx(end))-dy_s(idx(1)));
        duration_ps   = (numel(idx)-1) / max(1, fps_eff);
        drift_rate_ps = net_drift_ps / max(1e-9, duration_ps);
    else
        [med_step_ps,p95_step_ps,cum_path_ps,net_drift_ps,drift_rate_ps] = deal(NaN);
    end

    summaryRow = { ...
        files(f).name, nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, NaN, NaN, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, NaN, NaN, ...
        fps_eff, pxSize_um, eff_maxShift_um, defs.settle_min};
    
    if numel(summaryRow) ~= 20
        error('Summary row for %s has %d columns (expected 20).', files(f).name, numel(summaryRow));
    end
    allSummary(end+1,:) = summaryRow; %#ok<AGROW>
end

if ~isempty(allSummary)
    expectedCols = 20;
    if size(allSummary, 2) ~= expectedCols
        error('Expected %d columns in summary but got %d. Check compute_summary_from_csv.', ...
              expectedCols, size(allSummary, 2));
    end
    
    summaryTable = cell2table(allSummary, ...
        'VariableNames', {'File','n_frames', ...
        'MedianDisp_um','Pct95Disp_um','CumPath_um','NetDrift_um','DriftRate_um_per_s','Pegged_frac','Flatline_flag', ...
        'MedianDisp_postSettle_um','Pct95Disp_postSettle_um','CumPath_postSettle_um','NetDrift_postSettle_um','DriftRate_postSettle_um_per_s','Pegged_postSettle_frac','Flatline_postSettle_flag', ...
        'fps_used','px_um','max_shift_um','settle_min'});
    
    summaryTable.MATLAB_release   = repmat(string(versions_out.MATLAB_release), height(summaryTable), 1);
    summaryTable.NoRMCorre_commit = repmat(string(versions_out.NoRMCorre_commit), height(summaryTable), 1);
    
    writetable(summaryTable, fullfile(dataDir,'motion_summary.csv'));
    logf(logFid, 'Wrote summary: %s\n', fullfile(dataDir,'motion_summary.csv'));
end
end

% =====================================================================
% ==================== QC + ENRICHED (NEURO) ==========================
% =====================================================================
function qc_and_enriched_motion_folder_neuro(dataDir, cfg, defs_in, versions_in, logFid)
% ═══════════════════════════════════════════════════════════════════════════
% DUAL QC SYSTEM FOR NEURO DATA
% ═══════════════════════════════════════════════════════════════════════════
%
% This function produces TWO independent QC flags:
%
% ┌─────────────────────────────────────────────────────────────────────────┐
% │ MOVEMENT_VALID                                                          │
% │ ──────────────                                                          │
% │ Purpose: Can we trust motion measurements for immobilization inference? │
% │ Use for: Comparing Control vs Experimental movement metrics             │
% │ Checks:  Tracking integrity only (nan_frac, pegged, flatline, unstable) │
% │ Note:    Does NOT check boundary excursion (motion data still valid     │
% │          even if brain partially exits FOV)                             │
% └─────────────────────────────────────────────────────────────────────────┘
%
% ┌─────────────────────────────────────────────────────────────────────────┐
% │ ROI_ELIGIBLE                                                            │
% │ ────────────                                                            │
% │ Purpose: Is this file suitable for ROI detection and trace extraction?  │
% │ Use for: Deciding which files to run through EZcalcium/neural pipeline  │
% │ Checks:  MOVEMENT_VALID criteria PLUS boundary excursion check          │
% │ Note:    Stricter than MOVEMENT_VALID (brain must stay in FOV)          │
% └─────────────────────────────────────────────────────────────────────────┘
%
% A file can have MOVEMENT_VALID=true but ROI_ELIGIBLE=false if the brain
% partially exited the FOV. The motion tracking is still valid for
% immobilization comparison, but neural data would be incomplete.

px_um_default = cfg.pixel_size_um;
fps_default   = cfg.fps;

base_settle_min = cfg.settle_base_min;
min_settle_min  = cfg.settle_min_min;
frac_settle     = cfg.settle_frac;

% ═══════════════════════════════════════════════════════════════════════════
% MOVEMENT_VALID THRESHOLDS (tracking integrity for immobilization inference)
% ═══════════════════════════════════════════════════════════════════════════
max_pegged_frac       = cfg.qc_max_pegged_frac;
max_nan_frac          = cfg.qc_nan_frac_max;
% flatline and unstable_tracking are always fail conditions

% ═══════════════════════════════════════════════════════════════════════════
% ROI_ELIGIBLE THRESHOLDS (additional checks for neural pipeline)
% ═══════════════════════════════════════════════════════════════════════════
% boundary_fail_frac: If total excursion exceeds this fraction of FOV,
%                     brain may have exited and ROIs would be invalid

% Informational thresholds (not used for either QC, just for reporting)
max_spike_frac_info   = cfg.qc_max_spike_frac;
good_post_med_um_info = cfg.qc_good_post_med_um;

fps_override     = cfg.fps_override;
useUnicodeLabels = cfg.use_unicode_labels;
spike_thr_rate_um_s = cfg.qc_spike_thr_um_s;
good_post_med_um = cfg.qc_good_post_med_um;
spike_min_consec_user = cfg.qc_spike_min_consec_fr;
reposition_thr_um = cfg.reposition_thr_um;
reposition_min_sustain_s = cfg.reposition_min_sustain_s;
reposition_baseline_shift_um = cfg.reposition_baseline_shift_um;
boundary_warning_frac = cfg.boundary_warning_frac;
boundary_fail_frac = cfg.boundary_fail_frac;

umtxt = get_um_label(useUnicodeLabels);

% Derive CSV glob from file_pattern: replace extension with _shifts.csv
[~, pat_name, ~] = fileparts(cfg.file_pattern);
csv_pattern = [pat_name '_shifts.csv'];
csvs = dir(fullfile(dataDir, csv_pattern));
if isempty(csvs)
    logf(logFid, '[WARN] QC: no per-movie CSVs found in %s. Nothing to do.\n', dataDir);
    return
end

logf(logFid, '\n─── QC PROCESSING ───\n');
logf(logFid, 'Dual QC system: MOVEMENT_VALID (immobilization) + ROI_ELIGIBLE (neural pipeline)\n\n');

rows_basic = {}; rows_qc = {}; rows_enr = {};
tc_all = {}; psd_all = {};

for i = 1:numel(csvs)
    inCSV = fullfile(csvs(i).folder, csvs(i).name);
    [~,stem] = fileparts(inCSV);
    stem = regexprep(stem, '_shifts$', '');
    T = readtable(inCSV);

    px_um = px_um_default; if ismember('px_um', T.Properties.VariableNames), px_um = T.px_um(1); end
    
    % Determine fps
    if ~isempty(fps_override) && fps_override > 0
        fps = fps_override;
        fps_from_csv = fps_default;
        if ismember('fps_measured', T.Properties.VariableNames)
            fps_from_csv = T.fps_measured(1);
        elseif ismember('fps_used', T.Properties.VariableNames)
            fps_from_csv = T.fps_used(1);
        end
        if i == 1
            logf(logFid, '  Using fps_override=%.4f Hz (CSV reports %.1f Hz)\n', fps, fps_from_csv);
        end
    else
        fps = fps_default;   
        if ismember('fps_measured', T.Properties.VariableNames)
            fps = T.fps_measured(1);
        elseif ismember('fps_used', T.Properties.VariableNames)
            fps = T.fps_used(1);
        end
    end
    
    ny = fps/2;
    band_low  = [0.01, min(0.05, 0.2*ny)];
    band_mid  = [min(0.05, 0.2*ny), min(0.2, 0.8*ny)];
    band_high = [min(0.2, 0.8*ny), ny];
    
    spike_thr_um = spike_thr_rate_um_s / max(1,fps);
    if isempty(spike_min_consec_user)
        spike_min_consec_eff = max(6, round(0.18*fps));
    else
        spike_min_consec_eff = spike_min_consec_user;
    end

    aviPath = fullfile(csvs(i).folder, [stem '.avi']);
    if exist(aviPath,'file')==2
        vr = VideoReader(aviPath); H = vr.Height; W = vr.Width;
    else, H = 1024; W = 1024;
    end
    
    rescue_count = 0;
    unstable_tracking = false;
    if ismember('rescue_count', T.Properties.VariableNames)
        rescue_count = T.rescue_count(1);
    end
    if ismember('unstable_tracking', T.Properties.VariableNames)
        unstable_tracking = logical(T.unstable_tracking(1));
    end
    
    max_shift_px = compute_max_shift_px(H, W, px_um, cfg.max_shift_mode, cfg.max_shift_value);
    if isscalar(max_shift_px)
        maxShift_um = max_shift_px * px_um;
    else
        maxShift_um = max(max_shift_px) * px_um;
    end

    if all(ismember({'dx_um','dy_um'}, T.Properties.VariableNames))
        dx_um = T.dx_um; dy_um = T.dy_um;
    else
        dx_um = T.dx_px * px_um; dy_um = T.dy_px * px_um;
    end
    nF = numel(dx_um);
    steps_um = hypot([0; diff(dx_um)], [0; diff(dy_um)]);
    
    flat_eps_um = max(0.05, 0.25*px_um);
    
    % ════════════════════════════════════════════════════════════════════
    % BOUNDARY EXCURSION CHECK (for ROI_ELIGIBLE)
    % ════════════════════════════════════════════════════════════════════
    % Track total excursion to detect if brain may have exited FOV
    FOV_um = min(H, W) * px_um;  % ~266 µm for your setup
    
    dx_range_um = max(dx_um) - min(dx_um);  % Total X excursion
    dy_range_um = max(dy_um) - min(dy_um);  % Total Y excursion
    excursion_max_um = max(dx_range_um, dy_range_um);
    excursion_frac_FOV = excursion_max_um / FOV_um;
    
    boundary_warning = excursion_frac_FOV > boundary_warning_frac;
    boundary_exit = excursion_frac_FOV > boundary_fail_frac;
    
    if boundary_warning
        logf(logFid, '  [%s] Boundary warning: excursion %.1f%% of FOV (%.1f %s)\n', ...
            stem, excursion_frac_FOV*100, excursion_max_um, umtxt);
    end
    if boundary_exit
        logf(logFid, '  [%s] BOUNDARY EXIT: excursion %.1f%% of FOV - ROI_ELIGIBLE=false\n', ...
            stem, excursion_frac_FOV*100);
    end
    
    % Detect repositioning events
    reposition_mask = detect_reposition_events(steps_um, dx_um, dy_um, fps, ...
        reposition_thr_um, reposition_min_sustain_s, reposition_baseline_shift_um);
    
    if any(reposition_mask)
        reposMaskTable = table((1:nF)', double(reposition_mask), ...
            'VariableNames', {'Frame', 'IsReposition'});
        writetable(reposMaskTable, fullfile(csvs(i).folder, [stem '_reposition_mask.csv']));
        logf(logFid, '  [%s] Wrote reposition mask (%d frames)\n', stem, sum(reposition_mask));
    end
    
    thr_x_um = NaN; thr_y_um = NaN;
    if all(ismember({'max_shift_um_y','max_shift_um_x'}, T.Properties.VariableNames))
        thr_y_um = T.max_shift_um_y(1);
        thr_x_um = T.max_shift_um_x(1);
    elseif ismember('max_shift_um', T.Properties.VariableNames)
        thr_y_um = T.max_shift_um(1);
        thr_x_um = T.max_shift_um(1);
    else
        ms_px = compute_max_shift_px(H, W, px_um, cfg.max_shift_mode, cfg.max_shift_value);
        if isscalar(ms_px)
            thr_y_um = ms_px*px_um; 
            thr_x_um = thr_y_um;
        else
            thr_y_um = ms_px(1)*px_um; 
            thr_x_um = ms_px(2)*px_um;
        end
    end
    is_pegged = (abs(dx_um) > 0.98*thr_x_um) | (abs(dy_um) > 0.98*thr_y_um);
    pegged    = mean(is_pegged);
    nan_frac  = mean(~isfinite(dx_um) | ~isfinite(dy_um));
    
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
    duration_s = (nF - 1) / max(1, fps);
    duration_min = duration_s / 60;
    drift_rate = net_drift / max(1e-9, duration_s);

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

    settle_min_eff = min(base_settle_min, max(min_settle_min, frac_settle*duration_min));
    settle_frames  = min(nF, round(settle_min_eff * 60 * fps));
    pre_idx  = 1:settle_frames;
    post_idx = (settle_frames+1):nF;
    
    if exist('reposition_mask', 'var') && any(reposition_mask)
        post_idx_clean = post_idx(~reposition_mask(post_idx));
    else
        post_idx_clean = post_idx;
    end

    med_step_pre = median(steps_um(pre_idx),'omitnan');

    if numel(post_idx_clean) >= 2
        steps_post   = steps_um(post_idx_clean);
        med_step_ps  = median(steps_post,'omitnan');
        mean_step_ps = mean(steps_post,'omitnan');
        p95_step_ps  = Pctl(steps_post,95);
        p99_step_ps  = Pctl(steps_post,99);
        cum_path_ps  = sum(steps_post,'omitnan');

        net_drift_ps  = hypot(dx_s(post_idx_clean(end))-dx_s(post_idx_clean(1)), ...
                              dy_s(post_idx_clean(end))-dy_s(post_idx_clean(1)));
        dur_ps        = (numel(post_idx_clean)-1)/max(1,fps);
        drift_rate_ps = net_drift_ps / max(1e-9, dur_ps);

        pegged_ps   = mean(abs(dx_um(post_idx_clean)) > 0.98*thr_x_um | ...
                          abs(dy_um(post_idx_clean)) > 0.98*thr_y_um);
        flatline_ps = (std(dx_um(post_idx_clean)) < flat_eps_um) && ...
                      (std(dy_um(post_idx_clean)) < flat_eps_um);

        spike_frac_ps = frac_above_threshold(steps_post, spike_thr_um, spike_min_consec_eff);

        jitter_ratio = cum_path_ps / max(1e-9, net_drift_ps);
        jitter_ratio_display = min(jitter_ratio, 1000);

        [f_psd, P] = compute_psd_safe(steps_post, fps);
        [dom_f, dom_pow]   = dominant_peak(f_psd, P, [0.01 0.50]);
        if dom_f > 0
            dom_period_min = (1/dom_f) / 60;
        else
            dom_period_min = NaN;
        end
        frac_low  = band_power_frac(f_psd, P, band_low);
        frac_mid  = band_power_frac(f_psd, P, band_mid);
        frac_high = band_power_frac(f_psd, P, band_high);
        psd_all{end+1} = struct('f',f_psd,'P',P); %#ok<AGROW>
    else
        [med_step_ps,mean_step_ps,p95_step_ps,p99_step_ps,cum_path_ps,net_drift_ps,drift_rate_ps,pegged_ps,flatline_ps,spike_frac_ps,jitter_ratio,jitter_ratio_display] = deal(NaN);
        [dom_f, dom_pow, dom_period_min, frac_low, frac_mid, frac_high] = deal(NaN);
    end

    pre_floor_um = 0.02;
    if (duration_min < 2*settle_min_eff) || ~isfinite(med_step_ps) || ~isfinite(med_step_pre) || (med_step_pre < pre_floor_um)
        settle_ratio = NaN;
    else
        settle_ratio = med_step_ps / med_step_pre;
    end

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

    % ════════════════════════════════════════════════════════════════════
    % MOVEMENT_VALID: Tracking integrity for immobilization inference
    % ════════════════════════════════════════════════════════════════════
    % Can we trust these motion measurements for comparing immobilization
    % techniques? Only fails if tracking itself is unreliable.
    
    movement_fail_reasons = strings(0,1);
    
    if nan_frac >= max_nan_frac
        movement_fail_reasons(end+1) = "nan_frac_high";
    end
    
    if pegged >= max_pegged_frac
        movement_fail_reasons(end+1) = "tracking_ceiling";
    end
    
    if unstable_tracking
        movement_fail_reasons(end+1) = "unstable_tracking";
    end
    
    if flatline
        movement_fail_reasons(end+1) = "flatline_suspicious";
    end
    
    movement_valid = isempty(movement_fail_reasons);
    
    if isempty(movement_fail_reasons)
        movement_qc_reason = "TRACKING_OK";
    else
        movement_qc_reason = strjoin(movement_fail_reasons, ",");
    end
    
    % ════════════════════════════════════════════════════════════════════
    % ROI_ELIGIBLE: Suitable for neural ROI detection pipeline
    % ════════════════════════════════════════════════════════════════════
    % Requires MOVEMENT_VALID plus brain stayed in FOV
    
    roi_fail_reasons = movement_fail_reasons;  % Start with all movement failures
    
    if boundary_exit
        roi_fail_reasons(end+1) = "boundary_exit";
    end
    
    roi_eligible = isempty(roi_fail_reasons);
    
    if isempty(roi_fail_reasons)
        roi_qc_reason = "ELIGIBLE";
    else
        roi_qc_reason = strjoin(roi_fail_reasons, ",");
    end
    
    % ════════════════════════════════════════════════════════════════════
    % Health flags for quick visual inspection
    % ════════════════════════════════════════════════════════════════════
    
    % Movement health (for immobilization inference)
    MovementHealthFlag = "OK";
    if unstable_tracking
        MovementHealthFlag = "UNSTABLE_TRACKING";
    elseif flatline
        MovementHealthFlag = "FLATLINE";
    elseif ~movement_valid
        if nan_frac >= max_nan_frac, MovementHealthFlag = "NAN_EXCESS";
        elseif pegged >= max_pegged_frac, MovementHealthFlag = "TRACKING_CEILING";
        end
    end
    
    % ROI health (for neural pipeline)
    ROIHealthFlag = MovementHealthFlag;  % Start with movement health
    if movement_valid && boundary_exit
        ROIHealthFlag = "BOUNDARY_EXIT";  % Override if only boundary is the issue
    end
    
    % Motion quality (informational only)
    motion_quality = "GOOD";
    if isfinite(med_step_ps) && med_step_ps >= good_post_med_um_info
        motion_quality = "HIGH_MOTION";
    elseif isfinite(spike_frac_ps) && spike_frac_ps >= max_spike_frac_info
        motion_quality = "SPIKY";
    end

    % ════════════════════════════════════════════════════════════════════
    % Build output rows
    % ════════════════════════════════════════════════════════════════════

    rows_basic(end+1,:) = { ...
        [stem '.avi'], nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, pegged, flatline, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, pegged_ps, flatline_ps, ...
        fps, px_um, maxShift_um, settle_min_eff }; %#ok<AGROW>

    rows_qc(end+1,:) = { ...
        [stem '.avi'], nF, ...
        med_step, p95_step, cum_path, net_drift, drift_rate, pegged, flatline, ...
        med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, pegged_ps, flatline_ps, ...
        fps, px_um, maxShift_um, settle_min_eff, ...
        med_step_pre, settle_ratio, spike_frac_ps, jitter_ratio_display, ...
        nan_frac, rescue_count, logical(unstable_tracking), ...
        excursion_max_um, excursion_frac_FOV, logical(boundary_warning), logical(boundary_exit), ...
        logical(movement_valid), string(movement_qc_reason), string(MovementHealthFlag), ...
        logical(roi_eligible), string(roi_qc_reason), string(ROIHealthFlag), ...
        string(motion_quality)}; %#ok<AGROW>

    rows_enr(end+1,:) = { ...
        [stem '.avi'], nF, ...
        mean_step, sd_step, p99_step, ...
        mean_step_ps, p99_step_ps, ...
        pct_below_0p10, pct_below_0p20, ...
        drift_angle_deg, auto_settle_min, ...
        burst_count(steps_post_if(post_idx_clean,steps_um), 0.2, fps), burst_rate(steps_post_if(post_idx_clean,steps_um), 0.2, fps), ...
        burst_count(steps_post_if(post_idx_clean,steps_um), 0.5, fps), burst_rate(steps_post_if(post_idx_clean,steps_um), 0.5, fps), ...
        burst_count(steps_post_if(post_idx_clean,steps_um), 1.0, fps), burst_rate(steps_post_if(post_idx_clean,steps_um), 1.0, fps), ...
        dom_f, dom_period_min, dom_pow, ...
        frac_low, frac_mid, frac_high, ...
        med_step, med_step_ps, cum_path, cum_path_ps, ...
        pegged, pegged_ps, spike_frac_ps, jitter_ratio_display, ...
        excursion_max_um, excursion_frac_FOV, ...
        logical(movement_valid), string(MovementHealthFlag), ...
        logical(roi_eligible), string(ROIHealthFlag), ...
        rescue_count, logical(unstable_tracking)}; %#ok<AGROW>

    [min_idx, per_min_median] = per_minute_median(steps_um, fps);
    
    if numel(post_idx_clean) >= 2
        steps_post_full = steps_um(post_idx_clean);
        [min_idx_post, per_min_median_post] = per_minute_median(steps_post_full, fps);
        tc_all{end+1} = struct('minute',min_idx,'file_median',per_min_median, ...
                              'minute_post',min_idx_post,'file_median_post',per_min_median_post); %#ok<AGROW>
    else
        tc_all{end+1} = struct('minute',min_idx,'file_median',per_min_median, ...
                              'minute_post',[],'file_median_post',[]); %#ok<AGROW>
    end
end

% Create dynamic spike header
thr_txt = regexprep(sprintf('%.3g', spike_thr_rate_um_s), '\.', 'p');
spike_hdr = ['SpikeFrac_postSettle_over_' thr_txt 'um_per_s'];

hdr_basic = {'File','n_frames', ...
    'MedianDisp_um','Pct95Disp_um','CumPath_um','NetDrift_um','DriftRate_um_per_s','Pegged_frac','Flatline_flag', ...
    'MedianDisp_postSettle_um','Pct95Disp_postSettle_um','CumPath_postSettle_um','NetDrift_postSettle_um','DriftRate_postSettle_um_per_s','Pegged_postSettle_frac','Flatline_postSettle_flag', ...
    'fps_used','px_um','max_shift_um','settle_min'};

hdr_qc = [hdr_basic, ...
    {'MedianDisp_preSettle_um','SettleRatio_post_over_pre', spike_hdr, 'JitterRatio_post_cumPath_over_netDrift', ...
     'nan_frac','rescue_count','unstable_tracking', ...
     'Excursion_max_um','Excursion_frac_FOV','Boundary_warning','Boundary_exit', ...
     'MOVEMENT_VALID','Movement_QC_reason','Movement_HealthFlag', ...
     'ROI_ELIGIBLE','ROI_QC_reason','ROI_HealthFlag', ...
     'MotionQuality'}];

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
    spike_hdr,'JitterRatio_post_cumPath_over_netDrift', ...
    'Excursion_max_um','Excursion_frac_FOV', ...
    'MOVEMENT_VALID','Movement_HealthFlag', ...
    'ROI_ELIGIBLE','ROI_HealthFlag', ...
    'rescue_count','unstable_tracking'};

perfile_basic = cell2table(rows_basic,'VariableNames',hdr_basic);
perfile_qc    = cell2table(rows_qc,   'VariableNames',hdr_qc);
perfile_enr   = cell2table(rows_enr,  'VariableNames',hdr_enr);

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

% ════════════════════════════════════════════════════════════════════
% Group aggregates with dual QC summary
% ════════════════════════════════════════════════════════════════════
G = table;
G.n_files             = height(perfile_qc);
G.n_movement_valid    = sum(perfile_qc.MOVEMENT_VALID);
G.movement_valid_rate = G.n_movement_valid / max(1,G.n_files);
G.n_roi_eligible      = sum(perfile_qc.ROI_ELIGIBLE);
G.roi_eligible_rate   = G.n_roi_eligible / max(1,G.n_files);
G.n_boundary_exit     = sum(perfile_qc.Boundary_exit);
m = @(x) mean(x,'omitnan'); s = @(x) std(x,'omitnan');
G.mean_MedianDisp_post_um   = m(perfile_qc.MedianDisp_postSettle_um);
G.sd_MedianDisp_post_um     = s(perfile_qc.MedianDisp_postSettle_um);
G.mean_Excursion_frac_FOV   = m(perfile_qc.Excursion_frac_FOV);
G.max_Excursion_frac_FOV    = max(perfile_qc.Excursion_frac_FOV);
G.mean_rescue_count         = m(perfile_qc.rescue_count);
G.frac_unstable             = mean(perfile_qc.unstable_tracking);
writetable(G, fullfile(dataDir,'motion_group_summary.csv'));
logf(logFid, 'QC/ENR: wrote motion_group_summary.csv\n');

% Log dual QC summary
logf(logFid, '\n═══ DUAL QC SUMMARY ═══\n');
logf(logFid, 'Total files: %d\n', G.n_files);
logf(logFid, '  MOVEMENT_VALID: %d (%.1f%%) - usable for immobilization comparison\n', ...
    G.n_movement_valid, G.movement_valid_rate*100);
logf(logFid, '  ROI_ELIGIBLE:   %d (%.1f%%) - usable for neural pipeline\n', ...
    G.n_roi_eligible, G.roi_eligible_rate*100);
if G.n_boundary_exit > 0
    logf(logFid, '  Boundary exits: %d (MOVEMENT_VALID but not ROI_ELIGIBLE)\n', G.n_boundary_exit);
end
logf(logFid, '═══════════════════════\n\n');

% Per-minute timecourse
[tc_table, tc_plot, tc_table_post] = build_group_timecourse(tc_all);
writetable(tc_table, fullfile(dataDir,'group_timecourse.csv'));
logf(logFid, 'QC/ENR: wrote group_timecourse.csv\n');

rep_duration_min = median(perfile_basic.n_frames ./ perfile_basic.fps_used / 60, 'omitnan');
rep_settle_min = median(perfile_basic.settle_min, 'omitnan');
if ~isempty(tc_table_post) && rep_duration_min >= 2*rep_settle_min
    writetable(tc_table_post, fullfile(dataDir,'group_timecourse_post.csv'));
    logf(logFid, 'QC/ENR: wrote group_timecourse_post.csv\n');
end

set_figure_defaults();

rep_fps = max(1, median(perfile_qc.fps_used,'omitnan'));
if ~isfinite(rep_fps), rep_fps = 33.0; end
if isempty(spike_min_consec_user)
    spike_min_consec_for_label = max(6, round(0.18*rep_fps));
else
    spike_min_consec_for_label = spike_min_consec_user;
end

% ════════════════════════════════════════════════════════════════════
% DUAL QC STATUS FIGURE
% ════════════════════════════════════════════════════════════════════
try
    f_dual = figure('Visible','off','Color','w','Position',[100 100 1200 600]);
    
    tbl = sortrows(perfile_qc, {'ROI_ELIGIBLE','MOVEMENT_VALID'}, {'descend','descend'});
    y = 1:height(tbl);
    
    subplot(1,2,1);
    c1 = double(tbl.MOVEMENT_VALID);
    scatter(double(tbl.MOVEMENT_VALID), y, 60, c1, 'filled'); hold on;
    colormap(gca, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    set(gca,'YDir','reverse','YTick',y,'YTickLabel',tbl.File,'XTick',[0 1],'XTickLabel',{'INVALID','VALID'});
    xlim([-0.5 1.5]); ylim([0 height(tbl)+1]); grid on; box on;
    title({'MOVEMENT\_VALID','(for immobilization comparison)'});
    
    subplot(1,2,2);
    c2 = double(tbl.ROI_ELIGIBLE);
    scatter(double(tbl.ROI_ELIGIBLE), y, 60, c2, 'filled'); hold on;
    colormap(gca, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    set(gca,'YDir','reverse','YTick',y,'YTickLabel',tbl.File,'XTick',[0 1],'XTickLabel',{'INELIGIBLE','ELIGIBLE'});
    xlim([-0.5 1.5]); ylim([0 height(tbl)+1]); grid on; box on;
    title({'ROI\_ELIGIBLE','(for neural pipeline)'});
    
    sgtitle('Dual QC Status');
    
    exportgraphics(f_dual, fullfile(dataDir,'qc_dual_status.png'), 'Resolution', 300);
    exportgraphics(f_dual, fullfile(dataDir,'qc_dual_status.pdf'));
    close(f_dual);
    logf(logFid, 'Wrote qc_dual_status.png/pdf\n');
catch ME
    logf(logFid, '[WARN] Could not create dual QC figure: %s\n', ME.message);
end

% ════════════════════════════════════════════════════════════════════
% BOUNDARY EXCURSION FIGURE
% ════════════════════════════════════════════════════════════════════
try
    f_bound = figure('Visible','off','Color','w','Position',[100 100 800 500]);
    
    subplot(2,1,1);
    histogram(perfile_qc.Excursion_frac_FOV * 100, 'BinWidth', 5);
    hold on;
    xline(boundary_warning_frac*100, 'y--', 'Warning', 'LineWidth', 1.5);
    xline(boundary_fail_frac*100, 'r--', 'Fail', 'LineWidth', 1.5);
    xlabel('FOV Excursion (%)'); ylabel('Count');
    title('Boundary Excursion Distribution');
    grid on;
    
    subplot(2,1,2);
    c = double(perfile_qc.ROI_ELIGIBLE);
    scatter(perfile_qc.Excursion_frac_FOV * 100, perfile_qc.MedianDisp_postSettle_um, 50, c, 'filled');
    colormap(gca, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    hold on;
    xline(boundary_fail_frac*100, 'r--', 'LineWidth', 1.5);
    xlabel('FOV Excursion (%)'); ylabel(['Median step (' umtxt ')']);
    title('Excursion vs Motion (color = ROI\_ELIGIBLE)');
    grid on;
    
    exportgraphics(f_bound, fullfile(dataDir,'qc_boundary_excursion.png'), 'Resolution', 300);
    exportgraphics(f_bound, fullfile(dataDir,'qc_boundary_excursion.pdf'));
    close(f_bound);
    logf(logFid, 'Wrote qc_boundary_excursion.png/pdf\n');
catch ME
    logf(logFid, '[WARN] Could not create boundary figure: %s\n', ME.message);
end

% Overview figs (QC) - PNG + PDF
try
    f = figure('Visible','off','Color','w','Position',[100 100 980 560]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    c = double(perfile_qc.MOVEMENT_VALID);
    scatter(perfile_qc.MedianDisp_preSettle_um, perfile_qc.MedianDisp_postSettle_um, 36, c, 'filled');
    colormap(gca, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    hold on; yline(good_post_med_um,'k:'); xline(good_post_med_um,'k:');
    xlabel(['pre-settle median step (' umtxt ')']); ylabel(['post-settle median step (' umtxt ')']);
    title('Pre vs Post (color = MOVEMENT\_VALID)'); grid on;

    nexttile;
    histogram(perfile_qc.SettleRatio_post_over_pre, 'BinWidth',0.05);
    xlabel('post/pre median step'); ylabel('count'); title('Settle ratio'); grid on;

    nexttile;
    histogram(perfile_qc.(spike_hdr), 'BinWidth',0.005);
    xline(max_spike_frac_info,'r:','LineWidth',0.5); 
    xlabel(sprintf('spike fraction (>%.1f %s/s, ≥%d fr)', ...
        spike_thr_rate_um_s, umtxt, spike_min_consec_for_label));
    ylabel('count'); title('Spikes (informational)'); grid on;

    nexttile;
    histogram(perfile_qc.Pegged_frac, 'BinWidth',0.005);
    xline(max_pegged_frac,'r:'); xlabel('pegged fraction'); ylabel('count');
    title('Near-ceiling frames'); grid on;

    exportgraphics(f, fullfile(dataDir,'qc_overview.png'), 'Resolution', 300);
    exportgraphics(f, fullfile(dataDir,'qc_overview.pdf'));
    close(f);
catch
end

% Movement spectrum overview
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
        title('Group step power spectrum (post-settle)');
        
        ny = rep_fps/2;
        band_low  = [0.01, min(0.05, 0.2*ny)];
        band_mid  = [band_low(2), min(0.2, 0.8*ny)];
        band_high = [band_mid(2), ny];
        
        xline(band_low(2),':','Low','LabelHorizontalAlignment','left');
        xline(band_mid(2),':','Mid','LabelHorizontalAlignment','left');
        
        exportgraphics(f4, fullfile(dataDir,'movement_spectrum_overview.png'), 'Resolution', 300);
        exportgraphics(f4, fullfile(dataDir,'movement_spectrum_overview.pdf'));
        close(f4);
    end
catch
end

% Movement overview montage
try
    f5 = figure('Visible','off','Color','w','Position',[100 100 1100 800]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    try
        boxchart(ones(height(perfile_qc),1), perfile_qc.MedianDisp_postSettle_um);
    catch
        boxplot(perfile_qc.MedianDisp_postSettle_um);
    end
    hold on; yline(good_post_med_um,'r:');
    ylabel(['post-settle median step (' umtxt ')']); title('Across movies'); grid on;

    nexttile;
    try
        th = perfile_enr.DriftAngle_deg .* (pi/180);
        polarhistogram(th, 12, 'Normalization','probability'); title('Net drift direction');
    catch
        histogram(perfile_enr.DriftAngle_deg, -180:30:180); title('Net drift direction (deg)');
    end

    nexttile; hold on;
    x = tc_plot.minute; y = tc_plot.mean; e = tc_plot.sem;
    if ~isempty(x)
        fill([x; flipud(x)], [y-e; flipud(y+e)], [0.9 0.9 0.9], 'EdgeColor','none');
        plot(x, y, 'k', 'LineWidth', 1.4);
        xlabel('minute'); ylabel(['per-minute median step (' umtxt ')']); title('Group timecourse'); grid on;
    end

    nexttile;
    c = double(perfile_qc.MOVEMENT_VALID);
    xvals = max(perfile_qc.JitterRatio_post_cumPath_over_netDrift, 1e-6);
    yvals = perfile_qc.DriftRate_postSettle_um_per_s;
    
    ymax = max(yvals,[],'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax<=0, ymax = 1; end
    
    scatter(xvals, yvals, 40, c,'filled');
    colormap(gca, [0.75 0 0; 0 0.55 0]); caxis([-0.5 1.5]);
    set(gca, 'XScale', 'log');
    xlim([1e-6, 1e3]); ylim([0, ymax]);
    xlabel('JitterRatio [log]'); ylabel(['Drift rate (' umtxt '/s)']);
    title('Jitter vs drift (color = MOVEMENT\_VALID)'); grid on;

    exportgraphics(f5, fullfile(dataDir,'movement_overview.png'), 'Resolution', 300);
    exportgraphics(f5, fullfile(dataDir,'movement_overview.pdf'));
    close(f5);
catch
end

% ════════════════════════════════════════════════════════════════════
% HANDOFF FILES FOR TIFF PIPELINE (uses ROI_ELIGIBLE)
% ════════════════════════════════════════════════════════════════════
try
    stems = regexprep(perfile_qc.File, '\.avi$', '', 'ignorecase');
    
    post_start_fr = round(perfile_qc.settle_min .* 60 .* perfile_qc.fps_used) + 1;
    post_start_fr = min(post_start_fr, perfile_qc.n_frames);
    
    group_label = detect_group_label(dataDir);
    Group = repmat(string(group_label), height(perfile_qc), 1);
    
    roi_candidates = table( ...
        Group, ...
        stems, ...
        perfile_qc.File, ...
        perfile_qc.MOVEMENT_VALID, ...
        perfile_qc.Movement_QC_reason, ...
        perfile_qc.ROI_ELIGIBLE, ...
        perfile_qc.ROI_QC_reason, ...
        perfile_qc.Excursion_frac_FOV, ...
        post_start_fr, ...
        perfile_qc.n_frames, ...
        perfile_qc.fps_used, ...
        perfile_qc.px_um, ...
        'VariableNames', {'Group','Stem','AVI', ...
            'MOVEMENT_VALID','Movement_QC_reason', ...
            'ROI_ELIGIBLE','ROI_QC_reason', ...
            'Excursion_frac_FOV', ...
            'PostSettleStartFrame','n_frames','fps_used','px_um'});
    
    writetable(roi_candidates, fullfile(dataDir,'roi_candidates.csv'));
    logf(logFid, 'Pipeline handoff: wrote roi_candidates.csv\n');
    
    % JSON version
    try
        cand_struct = table2struct(roi_candidates);
        fid = fopen(fullfile(dataDir,'roi_candidates.json'),'w');
        if fid > 0
            fwrite(fid, jsonencode(cand_struct), 'char');
            fclose(fid);
        end
    catch
    end
    
    % Movement valid list (for immobilization comparison)
    try
        movement_valid_stems = stems(logical(perfile_qc.MOVEMENT_VALID));
        fid = fopen(fullfile(dataDir,'movement_valid_list.txt'),'w');
        if fid > 0
            for ii=1:numel(movement_valid_stems)
                fprintf(fid,'%s\n', movement_valid_stems{ii});
            end
            fclose(fid);
            logf(logFid, 'Wrote movement_valid_list.txt (%d files for immobilization comparison)\n', numel(movement_valid_stems));
        end
    catch
    end
    
    % ROI eligible list (for neural pipeline)
    try
        roi_eligible_stems = stems(logical(perfile_qc.ROI_ELIGIBLE));
        fid = fopen(fullfile(dataDir,'roi_eligible_list.txt'),'w');
        if fid > 0
            for ii=1:numel(roi_eligible_stems)
                fprintf(fid,'%s\n', roi_eligible_stems{ii});
            end
            fclose(fid);
            logf(logFid, 'Wrote roi_eligible_list.txt (%d files for neural pipeline)\n', numel(roi_eligible_stems));
        end
    catch
    end
    
catch ME
    logf(logFid, '[WARN] Could not write pipeline handoff files: %s\n', ME.message);
end

% ════════════════════════════════════════════════════════════════════
% Write thresholds.tsv with dual QC documentation
% ════════════════════════════════════════════════════════════════════
try
    thresholds = {
        % Header
        'DUAL_QC_SYSTEM', 'true', 'Two independent QC flags produced';
        
        % MOVEMENT_VALID thresholds
        'MOVEMENT_VALID_purpose', 'Immobilization inference', 'Can motion data be used to compare techniques?';
        'max_pegged_frac', max_pegged_frac, 'fraction (MOVEMENT_VALID threshold)';
        'max_nan_frac', max_nan_frac, 'fraction (MOVEMENT_VALID threshold)';
        'flatline_fail', true, 'boolean (MOVEMENT_VALID always fails)';
        'unstable_fail', true, 'boolean (MOVEMENT_VALID always fails)';
        
        % ROI_ELIGIBLE thresholds
        'ROI_ELIGIBLE_purpose', 'Neural pipeline suitability', 'Can file be used for ROI detection?';
        'boundary_warning_frac', boundary_warning_frac, 'fraction of FOV (warning only)';
        'boundary_fail_frac', boundary_fail_frac, 'fraction of FOV (ROI_ELIGIBLE fails)';
        
        % Informational thresholds
        'INFO_max_spike_frac', max_spike_frac_info, 'fraction (informational only)';
        'INFO_good_post_med_um', good_post_med_um_info, 'µm (informational only)';
        'spike_thr_rate_um_s', spike_thr_rate_um_s, 'µm/s';
        'spike_min_consec_fr', spike_min_consec_for_label, 'frames';
        'reposition_thr_um', reposition_thr_um, 'µm';
        'reposition_min_sustain_s', reposition_min_sustain_s, 's';
        'base_settle_min', base_settle_min, 'min';
        'representative_fps', rep_fps, 'Hz';
    };
    thrTbl = cell2table(thresholds, 'VariableNames', {'Parameter', 'Value', 'Units_or_Description'});
    writetable(thrTbl, fullfile(dataDir, 'thresholds.tsv'), 'FileType', 'text', 'Delimiter', '\t');
    logf(logFid, 'Wrote thresholds.tsv\n');
catch
end

% Write column dictionary
try
    write_column_dictionary(dataDir, boundary_warning_frac, boundary_fail_frac);
    logf(logFid, 'Wrote column_dictionary_motion.csv\n');
catch
end

% Write params JSON
try
    params = struct();
    params.analyzer       = "neuro_33fps_dual_qc";
    params.created_utc    = string(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
    if isfield(versions_in,'MATLAB_release'),   params.matlab_release    = string(versions_in.MATLAB_release);   else, params.matlab_release = "n/a"; end
    if isfield(versions_in,'NoRMCorre_commit'), params.normcorre_commit  = string(versions_in.NoRMCorre_commit); else, params.normcorre_commit = "n/a"; end
    
    params.fps_nominal    = 1000/30;
    params.fps_override   = fps_override;
    params.fps_representative = rep_fps;
    params.px_um          = px_um_default;
    
    rng_state = rng;
    params.rng_type = rng_state.Type;
    params.rng_seed = rng_state.Seed;

    % Dual QC documentation
    params.DUAL_QC_SYSTEM = struct( ...
        'description', 'Two independent QC flags for different purposes', ...
        'MOVEMENT_VALID', struct( ...
            'purpose', 'Can motion data be trusted for immobilization comparison?', ...
            'use_case', 'Comparing Control vs Experimental movement metrics', ...
            'checks', 'nan_frac, pegged_frac, flatline, unstable_tracking', ...
            'boundary_check', false), ...
        'ROI_ELIGIBLE', struct( ...
            'purpose', 'Is file suitable for ROI detection pipeline?', ...
            'use_case', 'Deciding which files to process through EZcalcium', ...
            'checks', 'All MOVEMENT_VALID checks PLUS boundary excursion', ...
            'boundary_check', true));
    
    params.MOVEMENT_VALID_thresholds = struct( ...
        'max_pegged_frac', max_pegged_frac, ...
        'max_nan_frac', max_nan_frac, ...
        'flatline_always_fail', true, ...
        'unstable_always_fail', true);
    
    params.ROI_ELIGIBLE_thresholds = struct( ...
        'inherits_from', 'MOVEMENT_VALID', ...
        'boundary_warning_frac', boundary_warning_frac, ...
        'boundary_fail_frac', boundary_fail_frac, ...
        'note', 'boundary_frac = max(dx_range, dy_range) / FOV_um');
    
    params.Informational_thresholds = struct( ...
        'note', 'NOT used for QC decisions, only for reporting', ...
        'good_post_med_um', good_post_med_um, ...
        'max_spike_frac', max_spike_frac_info, ...
        'spike_thr_rate_um_s', spike_thr_rate_um_s, ...
        'spike_min_consec_frames', spike_min_consec_for_label);

    params.Settle_policy = struct( ...
        'base_settle_min', base_settle_min, ...
        'min_settle_min', min_settle_min, ...
        'frac_settle', frac_settle);

    if ~isempty(defs_in)
        ncdefs = struct( ...
            'bin_width', defs_in.bin_width_default, ...
            'init_batch', defs_in.init_batch_default, ...
            'us_fac', defs_in.us_fac_default, ...
            'chunk_frames', defs_in.chunk_frames_default);
    else
        tmp = neuro_defaults(cfg);
        ncdefs = struct('bin_width',tmp.bin_width_default,'init_batch',tmp.init_batch_default,'us_fac',tmp.us_fac_default,'chunk_frames',tmp.chunk_frames_default);
    end
    params.NoRMCorre_defaults = ncdefs;

    jsonStr = jsonencode(params);
    fid = fopen(fullfile(dataDir,'params_motion_33fps.json'),'w');
    if fid>0
        fwrite(fid, jsonStr, 'char'); fclose(fid);
    end
catch ME
    logf(logFid, '[WARN] Could not write params JSON: %s\n', ME.message);
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
           '  analyze_motion_and_QC_neuro(''/path/to/data'')\n' ...
           '  analyze_motion_and_QC_neuro({''/path/A'', ''/path/B''})\n\n' ...
           'See also: neuro_motion_config']);
end
end

function cfg = merge_with_defaults(user_cfg)
% Merge a user-provided config struct with defaults. Unknown fields
% trigger a warning. Missing fields are filled from defaults.
cfg = neuro_motion_config();
default_fields = fieldnames(cfg);
user_fields = fieldnames(user_cfg);

for i = 1:numel(user_fields)
    if ismember(user_fields{i}, default_fields)
        cfg.(user_fields{i}) = user_cfg.(user_fields{i});
    else
        warning('neuro_motion:unknownConfig', ...
            'Unknown config field ''%s'' — ignored. See neuro_motion_config() for valid fields.', ...
            user_fields{i});
    end
end
end

function logf(fid, varargin)
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

function write_column_dictionary(dataDir, boundary_warning_frac, boundary_fail_frac)
dict = {
    % motion_qc.csv - Dual QC columns
    'motion_qc.csv', 'MOVEMENT_VALID', 'boolean', 'True if motion tracking is reliable for immobilization comparison';
    'motion_qc.csv', 'Movement_QC_reason', '', 'Reason for MOVEMENT_VALID failure (or TRACKING_OK)';
    'motion_qc.csv', 'Movement_HealthFlag', '', 'Summary health for movement data (OK, UNSTABLE_TRACKING, etc)';
    'motion_qc.csv', 'ROI_ELIGIBLE', 'boolean', 'True if file is suitable for ROI detection pipeline';
    'motion_qc.csv', 'ROI_QC_reason', '', 'Reason for ROI_ELIGIBLE failure (or ELIGIBLE)';
    'motion_qc.csv', 'ROI_HealthFlag', '', 'Summary health for neural pipeline (includes boundary check)';
    'motion_qc.csv', 'Excursion_max_um', 'µm', 'Maximum excursion in X or Y across recording';
    'motion_qc.csv', 'Excursion_frac_FOV', 'fraction', 'Excursion as fraction of FOV (brain exit risk)';
    'motion_qc.csv', 'Boundary_warning', 'boolean', sprintf('True if excursion > %.0f%% FOV (warning)', boundary_warning_frac*100);
    'motion_qc.csv', 'Boundary_exit', 'boolean', sprintf('True if excursion > %.0f%% FOV (ROI fails)', boundary_fail_frac*100);
    
    % Basic columns
    'motion_basic.csv', 'File', '', 'Filename of source AVI';
    'motion_basic.csv', 'n_frames', 'frames', 'Total number of frames processed';
    'motion_basic.csv', 'MedianDisp_um', 'µm', 'Median per-frame displacement (full movie)';
    'motion_basic.csv', 'MedianDisp_postSettle_um', 'µm', 'Median displacement (post-settle)';
    'motion_basic.csv', 'Pegged_frac', 'fraction', 'Fraction of frames at registration ceiling';
    'motion_basic.csv', 'Flatline_flag', 'boolean', 'True if signal is flat (no movement detected)';
    'motion_basic.csv', 'fps_used', 'Hz', 'Frame rate used for analysis';
    'motion_basic.csv', 'px_um', 'µm/px', 'Pixel size in micrometers';
    
    % QC columns
    'motion_qc.csv', 'nan_frac', 'fraction', 'Fraction of NaN or non-finite frames';
    'motion_qc.csv', 'rescue_count', 'count', 'Number of registration window rescues';
    'motion_qc.csv', 'unstable_tracking', 'boolean', 'True if tracking hit rescue cap repeatedly';
    'motion_qc.csv', 'MotionQuality', '', 'Informational: GOOD, HIGH_MOTION, or SPIKY';
    
    % Handoff files
    'roi_candidates.csv', 'MOVEMENT_VALID', 'boolean', 'Use for immobilization comparison';
    'roi_candidates.csv', 'ROI_ELIGIBLE', 'boolean', 'Use for neural pipeline (check this for ROI detection)';
    'movement_valid_list.txt', '', '', 'File stems passing MOVEMENT_VALID (for immobilization analysis)';
    'roi_eligible_list.txt', '', '', 'File stems passing ROI_ELIGIBLE (for neural pipeline)';
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

function defs = neuro_defaults(cfg)
% Build internal defaults struct from user config.
defs = struct();
defs.px_um               = cfg.pixel_size_um;
defs.fps_nominal         = cfg.fps;
defs.fps                 = cfg.fps;
defs.bin_width_default   = cfg.bin_width;
defs.us_fac_default      = cfg.us_fac;
defs.init_batch_default  = cfg.init_batch;
defs.chunk_frames_default= cfg.chunk_frames;
defs.settle_min          = cfg.settle_base_min;
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

function g = detect_group_label(p)
pl = lower(char(string(p)));
if contains(pl, 'hydrogel +') || contains(pl, 'ether') || contains(pl, 'experimental')
    g = "experimental";
elseif contains(pl, 'hydrogel') || contains(pl, 'control')
    g = "control";
else
    g = "unknown";
end
end

function max_shift_px = compute_max_shift_px(H, W, px_um, mode, val)
minHW = min(H,W);
half_lim = floor(minHW/2) - 1;
switch lower(string(mode))
    case "auto"
        max_shift_px = max(1, min(half_lim, round(0.40 * minHW)));
    case "auto_conservative"
        max_shift_px = max(1, min(half_lim, round(0.10 * minHW)));
    case "auto_loose"
        max_shift_px = max(1, min(half_lim, round(0.45 * minHW)));
    case "auto_extreme"
        max_shift_px = max(1, min(half_lim, round(0.50 * minHW)));
    case "fixed_um"
        if isempty(val), error('max_shift_value (um) required for fixed_um'); end
        v = val(:)'; if numel(v)==1, v=[v v]; end
        ms = round(v./px_um);
        max_shift_px = [max(1, min(ms(1), floor(H/2)-1)), max(1, min(ms(2), floor(W/2)-1))];
        if numel(unique(max_shift_px)) == 1
            max_shift_px = max_shift_px(1);
        end
    case "fixed_px"
        if isempty(val), error('max_shift_value (px) required for fixed_px'); end
        v = val(:)'; if numel(v)==1, v=[v v]; end
        ms = round(v);
        max_shift_px = [max(1, min(ms(1), floor(H/2)-1)), max(1, min(ms(2), floor(W/2)-1))];
        if numel(unique(max_shift_px)) == 1
            max_shift_px = max_shift_px(1);
        end
    otherwise
        error('Unknown max_shift_mode: %s', mode);
end
end

function rescue_cap = compute_rescue_cap(H, W, initial)
if isscalar(initial)
    rescue_cap = max(1, min(floor(min(H,W)/2)-1, round(0.50 * min(H,W))));
else
    rescue_cap = [max(1, min(floor(H/2)-1, round(0.50 * H))), ...
                  max(1, min(floor(W/2)-1, round(0.50 * W)))];
end
if isscalar(initial) && isscalar(rescue_cap)
    rescue_cap = max(rescue_cap, initial);
elseif ~isscalar(initial) && ~isscalar(rescue_cap)
    rescue_cap = max(rescue_cap, initial);
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
        [Pxx,f] = pwelch(x, hamming(min(2048, max(128, floor(numel(x)/4)))), [], [], fps);
        return;
    catch
    end
end
N = 2^nextpow2(max(128, numel(x)));
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

function [tc_table, tc_plot, tc_table_post] = build_group_timecourse(tc_all)
maxM = 0; 
maxM_post = 0;
for ii=1:numel(tc_all)
    maxM = max(maxM, max(tc_all{ii}.minute));
    if isfield(tc_all{ii}, 'minute_post') && ~isempty(tc_all{ii}.minute_post)
        maxM_post = max(maxM_post, max(tc_all{ii}.minute_post));
    end
end

M = (1:maxM).';
vals = NaN(maxM, numel(tc_all));
for j=1:numel(tc_all)
    m = tc_all{j}.minute; 
    v = tc_all{j}.file_median; 
    vals(m,j) = v; 
end
mean_v = mean(vals, 2, 'omitnan'); 
n_v = sum(~isnan(vals),2);
sem_v  = std(vals, 0, 2, 'omitnan') ./ sqrt(max(1,n_v));

tc_table_post = [];

if maxM_post > 0
    M_post = (1:maxM_post).';
    vals_post = NaN(maxM_post, numel(tc_all));
    for j=1:numel(tc_all)
        if isfield(tc_all{j}, 'minute_post') && ~isempty(tc_all{j}.minute_post)
            m_post = tc_all{j}.minute_post; 
            v_post = tc_all{j}.file_median_post;
            vals_post(m_post,j) = v_post;
        end
    end
    mean_v_post = mean(vals_post, 2, 'omitnan');
    n_v_post = sum(~isnan(vals_post),2);
    sem_v_post = std(vals_post, 0, 2, 'omitnan') ./ sqrt(max(1,n_v_post));
    
    tc_table = table(M, n_v, mean_v, sem_v, 'VariableNames', ...
        {'minute','n_files_active','mean_perMinuteMedian_um','sem_perMinuteMedian_um'});
    
    tc_table_post = table(M_post, n_v_post, mean_v_post, sem_v_post, 'VariableNames', ...
        {'minute_post_settle','n_files_active_post','mean_perMinuteMedian_post_um','sem_perMinuteMedian_post_um'});
    
    tc_plot = struct('minute',M,'mean',mean_v,'sem',sem_v, ...
                     'minute_post',M_post,'mean_post',mean_v_post,'sem_post',sem_v_post);
else
    tc_table = table(M, n_v, mean_v, sem_v, 'VariableNames', ...
        {'minute','n_files_active','mean_perMinuteMedian_um','sem_perMinuteMedian_um'});
    tc_plot = struct('minute',M,'mean',mean_v,'sem',sem_v);
end
end

function [dx, dy] = extract_dx_dy(shifts)
% Extract dx (horizontal/X) and dy (vertical/Y) from NoRMCorre shifts
% IMPORTANT: NoRMCorre uses [row, col] = [Y, X] convention
%
% OUTPUT CONVENTION (MOTION direction):
%   dx > 0 = sample moved RIGHT
%   dy > 0 = sample moved UP
%
% NoRMCorre has asymmetric sign conventions:
%   - X output is inverted, so we NEGATE
%   - Y output is already motion direction, keep as-is
dx = []; dy = [];
if isstruct(shifts)
    if isfield(shifts,'x') && isfield(shifts,'y')
        dx = -double(shifts.x(:));  % NEGATE X (was inverted)
        dy = double(shifts.y(:));   % Keep Y as-is
        return;
    end
    if isfield(shifts,'shifts')
        try
            % NoRMCorre convention: shifts(1)=row=Y, shifts(2)=col=X
            tmp = arrayfun(@(s) s.shifts(:), shifts, 'UniformOutput', false);
            tmp = cat(2, tmp{:}); 
            dx = -double(tmp(2,:)).';  % col shift = X, NEGATE (inverted)
            dy = double(tmp(1,:)).';   % row shift = Y, keep as-is
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
    dx = -double(val(:,2));   % column 2 = X, NEGATE
    dy = double(val(:,1));    % column 1 = Y, as-is
elseif numel(sz)==2 && sz(1)==2
    dx = -double(val(2,:)).'; % row 2 = X, NEGATE
    dy = double(val(1,:)).';  % row 1 = Y, as-is
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

function S = compute_summary_from_csv(csvPath, pxSize_um, fps_default, max_shift_px_in, settle_min)
T = readtable(csvPath);
if ismember('px_um', T.Properties.VariableNames),       pxSize_um   = T.px_um(1);     end
if ismember('fps_measured', T.Properties.VariableNames)
    fps_default = T.fps_measured(1);
elseif ismember('fps_used', T.Properties.VariableNames)
    fps_default = T.fps_used(1);
end
if ismember('max_shift_um', T.Properties.VariableNames)
    maxShift_um = T.max_shift_um(1);
elseif ismember('max_shift_px', T.Properties.VariableNames)
    maxShift_um = T.max_shift_px(1) * pxSize_um;
else
    if ~isnan(max_shift_px_in), maxShift_um = max_shift_px_in * pxSize_um; else, maxShift_um = NaN; end
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

settle_frames = min(nF, round(settle_min * 60 * fps_default));
idx = (settle_frames+1):nF;
if numel(idx) >= 2
    med_step_ps   = median(steps_um(idx), 'omitnan');
    p95_step_ps   = Pctl(steps_um(idx), 95);
    cum_path_ps   = sum(steps_um(idx), 'omitnan');
    dx_s2 = dx_s; dy_s2 = dy_s;
    net_drift_ps  = hypot(dx_s2(idx(end))-dx_s2(idx(1)), dy_s2(idx(end))-dy_s2(idx(1)));
    duration_ps   = (numel(idx)-1) / max(1, fps_default);
    drift_rate_ps = net_drift_ps / max(1e-9, duration_ps);
else
    [med_step_ps,p95_step_ps,cum_path_ps,net_drift_ps,drift_rate_ps] = deal(NaN);
end

[~, base] = fileparts(csvPath);
fileStem = regexprep(base, '_shifts$', '');
S = { [fileStem '.avi'], nF, ...
    med_step, p95_step, cum_path, net_drift, drift_rate, NaN, NaN, ...
    med_step_ps, p95_step_ps, cum_path_ps, net_drift_ps, drift_rate_ps, NaN, NaN, ...
    fps_default, pxSize_um, maxShift_um, settle_min };
end

function v = steps_post_if(post_idx, steps_um)
if isempty(post_idx) || numel(post_idx)<2, v = []; else, v = steps_um(post_idx); end
end

function [starts, ends] = find_runs(x)
x = x(:) ~= 0; 
d = diff([false; x; false]);
starts = find(d==1); 
ends = find(d==-1)-1;
end

function mask = detect_reposition_events(steps_um, dx_um, dy_um, fps, thr_um, min_sustain_s, baseline_shift_um)
if nargin < 5, thr_um = 2.0; end
if nargin < 6, min_sustain_s = 2.0; end
if nargin < 7, baseline_shift_um = 2.0; end

mask = false(size(steps_um));
hi = steps_um > thr_um;
[st, en] = find_runs(hi);

for k = 1:numel(st)
    if (en(k) - st(k) + 1) >= round(min_sustain_s*fps)
        pre_idx = max(1, st(k)-round(2*fps)) : max(1, st(k)-1);
        post_idx = min(en(k)+1, numel(steps_um)) : min(en(k)+round(2*fps), numel(steps_um));
        
        if numel(pre_idx) > 10 && numel(post_idx) > 10
            dx_pre = median(dx_um(pre_idx), 'omitnan');
            dy_pre = median(dy_um(pre_idx), 'omitnan');
            dx_post = median(dx_um(post_idx), 'omitnan');
            dy_post = median(dy_um(post_idx), 'omitnan');
            
            baseline_shift = hypot(dx_post - dx_pre, dy_post - dy_pre);
            if baseline_shift > baseline_shift_um
                event_start = max(1, st(k) - round(0.5*fps));
                event_end = min(numel(steps_um), en(k) + round(0.5*fps));
                mask(event_start:event_end) = true;
            end
        end
    end
end
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