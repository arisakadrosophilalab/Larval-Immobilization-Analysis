# Multimodal Immobilization of Drosophila Larvae for In Vivo Calcium Imaging

Analysis code for the manuscript: "Multimodal immobilization of second-instar *Drosophila melanogaster* larvae using PF-127 hydrogel and diethyl ether for calcium imaging" (submitted to the *Journal of Neuroscience Methods*).

## About This Code

All analysis code in this repository is computationally equivalent to the code used to generate the results reported in the manuscript. The code has been refactored for portability: hardcoded lab paths have been replaced with folder pickers and configuration files, parameters are documented and adjustable, and help documentation has been added. **No computational logic, formulas, thresholds, or statistical methods were changed.** The original lab-internal code is available from the authors upon request.

## Repository Structure

```
Larval-Immobilization-Analysis/
│
├── README.md
├── LICENSE
│
├── motion-analysis/
│   ├── analyze-motion/
│   │   ├── analyze_motion_and_QC.m          # Motion tracking for GFP recordings (1 Hz)
│   │   └── motion_analysis_config.m         # Configuration
│   │
│   └── compare-motion-groups/
│       ├── compare_motion_groups.m           # Statistical group comparison of motion
│       └── compare_motion_config.m           # Configuration
│
└── calcium-imaging/
    ├── run-this-1st/
    │   ├── analyze_motion_and_QC_neuro.m     # Preprocessing: extracts shifts + QC for pipeline
    │   └── neuro_motion_config.m             # Configuration
    │
    ├── run-this-2nd/
    │   ├── COMPLETE_PIPELINE_RUN_ME_v4_5.m   # ROI detection + trace extraction
    │   └── pipeline_config.m                 # Configuration
    │
    ├── run-this-3rd/
    │   ├── Assess_ROI_Quality.m              # Signal quality assessment
    │   └── assess_quality_config.m           # Configuration
    │
    ├── run-this-4th/
    │   ├── ROI_REFINEMENT_v3_2.m             # Multi-criteria ROI quality control
    │   └── roi_refinement_config.m           # Configuration
    │
    └── run-this-5th/
        ├── COMPARE_NEURO_CONDITIONS_v3.m     # Statistical comparison of conditions
        └── compare_conditions_config.m       # Configuration
```

## Which Scripts Do I Need?

This repository supports two workflows. Use whichever applies to your experiment.

### Workflow A: Movement Analysis Only

Use this if you want to quantify how well an immobilization method works by tracking residual motion in wide-field fluorescence recordings. This is the workflow that produces the motion comparison results in the paper (Tables 2–4, Figures 2–5).

**When to use:** You have GFP-expressing larvae recorded at low frame rate (~1 Hz) and want to measure and statistically compare how much they move under different immobilization conditions.

**Step 1 — Track motion in each condition's video folder:**

```matlab
cfg = motion_analysis_config();
cfg.pixel_size_um = 266/1024;          % your pixel size
cfg.fps = 1;                            % your frame rate
cfg.file_pattern = 'S*.avi';            % your file naming convention
analyze_motion_and_QC('/path/to/control/videos', cfg);
analyze_motion_and_QC('/path/to/experimental/videos', cfg);
```

**Step 2 — Compare groups statistically:**

```matlab
cfg = compare_motion_config();
cfg.control_label = 'Hydrogel';
cfg.experimental_label = 'Diethyl Ether + Hydrogel';
compare_motion_groups('/path/to/control', '/path/to/experimental', '/path/to/output', cfg);

% To analyze only the first 30 minutes (as in the manuscript's exploratory analysis):
cfg.time_window_min = [0 30];
compare_motion_groups('/path/to/control', '/path/to/experimental', '/path/to/output_30min', cfg);
```

### Workflow B: Calcium Imaging Pipeline

Use this if you want to detect neurons, extract calcium traces, and compare neural activity between conditions.

**When to use:** You have GCaMP-expressing larvae recorded at high frame rate (~33 Hz) as BigTIFF files and you want ROIs, calcium traces, and statistical comparisons between conditions.

The calcium imaging workflow has 5 steps. Each step consumes the output of the previous one:

```
Step 0: analyze_motion_and_QC_neuro.m       → shift CSVs + QC flags (ROI_ELIGIBLE list)
          ↓
Step 1: COMPLETE_PIPELINE_RUN_ME_v4_5.m     → motion-corrected TIFFs, ROIs, ΔF/F traces
          ↓
Step 2: Assess_ROI_Quality.m                → signal quality flags, sample eligibility lists
          ↓
Step 3: ROI_REFINEMENT_v3_2.m              → refined ROI sets (false positives removed)
          ↓
Step 4: COMPARE_NEURO_CONDITIONS_v3.m      → statistical comparison between conditions
```

**Step 0 — Motion tracking and QC gating (REQUIRED before the pipeline)**

This is a prerequisite for the calcium pipeline. You must run `analyze_motion_and_QC_neuro.m` on the AVI versions of your neural recordings before running the pipeline. This script:

- Performs NoRMCorre rigid registration on the AVI videos to estimate frame-by-frame shifts
- Saves per-frame shift CSVs (`*_shifts.csv`) that the pipeline reads in Step 1 to apply motion correction to the original BigTIFF files
- Produces a dual QC system that determines which files enter the pipeline:
  - **MOVEMENT_VALID** — Is the motion tracking reliable? (registration didn't fail)
  - **ROI_ELIGIBLE** — Is this file suitable for ROI detection? (brain stayed in the field of view)
- Generates an `roi_eligible_list.txt` that the pipeline uses to know which files to process

The pipeline will not run without the shift CSVs from this step.

```matlab
cfg = neuro_motion_config();
cfg.pixel_size_um = 266/1024;
cfg.fps_override = 33.33;
analyze_motion_and_QC_neuro('/path/to/control/avi_videos', cfg);
analyze_motion_and_QC_neuro('/path/to/experimental/avi_videos', cfg);
```

**Why AVI and TIFF?** The neural recordings exist in two formats: compressed AVI files (smaller, used for motion estimation) and raw BigTIFF files (larger, full quality, used for ROI detection). Step 0 processes the AVIs to compute shifts; Step 1 applies those shifts to the TIFFs.

**Step 1 — ROI detection and trace extraction**

`COMPLETE_PIPELINE_RUN_ME_v4_5.m` reads the shift CSVs from Step 0 and processes only the files listed in `roi_eligible_list.txt`. For each file it: applies rigid motion correction to the BigTIFF using Bio-Formats, generates a brain mask, computes local correlation (Cn) and peak-to-noise ratio (PNR) maps, detects ROI seeds via adaptive thresholding, grows ROI regions using correlation-based expansion, extracts raw fluorescence traces, applies neuropil subtraction (coefficient = 0.7), computes ΔF/F with a sliding baseline, optionally applies global signal regression to remove shared neuropil signal, smooths traces with a Savitzky-Golay filter, and detects calcium events using GCaMP8f-optimized peak finding.

```matlab
config = pipeline_config();
config.pixelSize_um = 266/1024;
config.frameRate = 33.33;
config.conditions(1).name = 'Control';
config.conditions(1).tifDir = '/path/to/control/tif';
config.conditions(1).aviDir = '/path/to/control/avi';
COMPLETE_PIPELINE_RUN_ME_v4_5(config);
```

**Step 2 — Signal quality assessment**

`Assess_ROI_Quality.m` evaluates signal quality across all samples to identify artifacts and neuropil contamination. It categorizes each sample based on pairwise synchrony:

| Category | Pairwise Synchrony | Status |
|---|---|---|
| OK | < 0.15 | Clean signal |
| MILD_NEUROPIL | 0.15 – 0.30 | Acceptable |
| HIGH_NEUROPIL | 0.30 – 0.50 | Needs regression |
| SEVERE_NEUROPIL | > 0.50 | Unusable |
| ARTIFACT | — | Extreme ΔF/F values (> 5 or < −2) |

```matlab
cfg = assess_quality_config();
% Edit condition paths in cfg, then:
Assess_ROI_Quality(cfg);
```

**Step 3 — ROI refinement**

`ROI_REFINEMENT_v3_2.m` applies multi-criteria quality control to remove false-positive ROIs. Criteria include:

- **Shape:** Edge contact (>20% perimeter on brain mask edge), compactness, bounding box fill (<35% or >95%), aspect ratio (>3), solidity (<0.5), eccentricity (>0.95)
- **Trace quality:** Negative transients, rise/decay kinetics inconsistent with GCaMP8f, single-frame step artifacts, sustained plateaus without decay
- **Statistics:** Pairwise trace correlation >0.7 (removes the lower-SNR duplicate), extreme SNR (>150), zero detected events, baseline drift

Each removed ROI is logged with the specific reason for rejection.

```matlab
cfg = roi_refinement_config();
% Edit condition paths in cfg, then:
ROI_REFINEMENT_v3_2(cfg);
```

**Step 4 — Statistical comparison of conditions**

`COMPARE_NEURO_CONDITIONS_v3.m` performs a publication-grade comparison of the refined calcium data (post-refinement, not raw pipeline output). Metrics compared include ROI detection quality (yield, removal rate, 11-category removal breakdown), signal quality (SNR, event rates, event amplitudes), motion-induced noise (pairwise synchrony, correlated ROI removal rate), and temporal stability (baseline drift). Statistics include Mann-Whitney U tests, bootstrapped 95% confidence intervals, Cliff's delta and Cohen's d effect sizes, and Bonferroni correction for multiple comparisons.

```matlab
cfg = compare_conditions_config();
% Edit condition paths to point to your Refined Data folders, then:
COMPARE_NEURO_CONDITIONS_v3(cfg);
```

## Why Are There Two Motion Scripts?

These two scripts share the same core NoRMCorre registration engine but serve completely different purposes:

**`analyze_motion_and_QC.m`** is the motion analysis script. It processes GFP movement recordings (1 Hz, ~60 min) and its outputs — motion metrics like mean speed, median step, and jitter ratio — are the actual data that get statistically compared between conditions using `compare_motion_groups.m`. This is what produces the movement results in the paper.

**`analyze_motion_and_QC_neuro.m`** is a preprocessing step for the calcium pipeline. It processes the AVI versions of neural recordings (33.33 Hz, ~5 min) to extract two things the pipeline needs: (1) per-frame shift CSVs that tell the pipeline how to motion-correct the BigTIFF files, and (2) QC flags (MOVEMENT_VALID, ROI_ELIGIBLE) that tell the pipeline which files are safe to process. Its motion data is never compared between groups — it exists only to set up the calcium pipeline.

| | analyze_motion_and_QC.m | analyze_motion_and_QC_neuro.m |
|---|---|---|
| **Role** | Analysis endpoint | Preprocessing for calcium pipeline |
| **Input** | GFP movement recordings (1 Hz) | GCaMP neural recordings, AVI versions (33.33 Hz) |
| **What happens to its output** | Compared between groups | Fed into COMPLETE_PIPELINE_RUN_ME_v4_5.m |
| **QC system** | Single flag (QC_PASS) | Dual flags (MOVEMENT_VALID + ROI_ELIGIBLE) |
| **Key outputs** | Motion metrics for statistical comparison | Shift CSVs + roi_eligible_list.txt for pipeline |

## Installation and Dependencies

### MATLAB

MATLAB R2020b or later is required (R2024b recommended). The following toolboxes are needed:

| Toolbox | Required by | Notes |
|---|---|---|
| Image Processing Toolbox | All scripts | imtranslate, regionprops, morphological operations |
| Statistics and Machine Learning Toolbox | Comparison scripts | ranksum, lillietest, prctile |
| Signal Processing Toolbox | Motion scripts | medfilt1, pwelch (FFT fallback available if missing) |

### NoRMCorre

NoRMCorre performs rigid motion registration. Required for both workflows.

```bash
git clone https://github.com/flatironinstitute/NoRMCorre.git
```

Then in MATLAB:

```matlab
addpath(genpath('/path/to/NoRMCorre'));
savepath;  % optional: save so you don't have to add it every session
```

### Bio-Formats (Calcium Pipeline Only)

Bio-Formats reads BigTIFF files in the calcium imaging pipeline. The motion analysis scripts do not need Bio-Formats.

1. Download the MATLAB toolbox from: https://www.openmicroscopy.org/bio-formats/downloads/
2. Unzip to a permanent location (e.g., `C:\toolboxes\bfmatlab\`)
3. Add to MATLAB path:
   ```matlab
   addpath(genpath('/path/to/bfmatlab'));
   savepath;
   ```
4. Verify:
   ```matlab
   reader = bfGetReader();  % should not error
   ```

### BigTIFF Writing

The calcium pipeline writes motion-corrected output as BigTIFF files (>4 GB). This uses MATLAB's built-in `Tiff` class with the `'w8'` flag (available since R2014b). No additional installation needed.

### Quick Dependency Check

```matlab
% NoRMCorre
assert(exist('normcorre', 'file') == 2, 'NoRMCorre not found. See README for installation.');

% Bio-Formats (only needed for calcium pipeline)
assert(exist('bfGetReader', 'file') == 2, 'Bio-Formats not found. See README for installation.');

% MATLAB Toolboxes
assert(exist('imtranslate', 'file') == 2, 'Image Processing Toolbox not found.');
assert(exist('prctile', 'file') == 2, 'Statistics and Machine Learning Toolbox not found.');
assert(exist('medfilt1', 'file') == 2, 'Signal Processing Toolbox not found.');

disp('All dependencies OK');
```

## Configuration

Each script accepts an optional configuration struct from its companion config file. Default parameters match those used in the manuscript. To customize:

```matlab
cfg = motion_analysis_config();    % load defaults
cfg.pixel_size_um = 0.5;           % change what you need
cfg.fps = 2;
analyze_motion_and_QC('/data', cfg);
```

If you omit the config argument, defaults are used and a folder picker dialog opens for path selection.

The calcium pipeline scripts each have a dedicated config file:

| Script | Config file |
|---|---|
| analyze_motion_and_QC_neuro.m | neuro_motion_config.m |
| COMPLETE_PIPELINE_RUN_ME_v4_5.m | pipeline_config.m |
| Assess_ROI_Quality.m | assess_quality_config.m |
| ROI_REFINEMENT_v3_2.m | roi_refinement_config.m |
| COMPARE_NEURO_CONDITIONS_v3.m | compare_conditions_config.m |

### Key Parameters to Set for Your System

| Parameter | Where | What it is |
|---|---|---|
| `pixel_size_um` | All config files | Pixel size in µm (sensor pixel pitch ÷ magnification) |
| `fps` | motion_analysis_config.m | Frame rate for GFP movement recordings |
| `fps_override` | neuro_motion_config.m | True acquisition rate for neural recordings |
| `file_pattern` | Config files | Glob pattern matching your video filenames (e.g., `'S*.avi'`) |
| `config.frameRate` | pipeline_config.m | Frame rate for calcium imaging (should match fps_override) |
| `config.pixelSize_um` | pipeline_config.m | Same pixel size as above |
| `config.conditions` | pipeline_config.m | Paths to your TIF and AVI data folders |

## Shift Convention

All scripts use a consistent motion-direction convention:

- `dx_px > 0` = sample moved right
- `dy_px > 0` = sample moved down

For correction, the pipeline applies `imtranslate(img, [-dx, +dy])`. The sign asymmetry between x and y is intentional: it arises from the conversion between Cartesian coordinates (y-up) used by the motion tracker and image coordinates (y-down) used by `imtranslate`. The column dictionary generated by the pipeline includes a detailed explanation.

## Output Files

### Motion Analysis Outputs (per folder)

| File | Description |
|---|---|
| `*_shifts.csv` | Per-frame shift estimates (dx, dy in pixels and µm) |
| `*_motion_trace.png/pdf` | Motion trace visualization |
| `motion_qc.csv` | QC flags and motion metrics for all files |
| `motion_enriched.csv` | Extended metrics (PSD, bursts, drift) |
| `motion_group_summary.csv` | Group-level aggregates |
| `roi_eligible_list.txt` | Files passing ROI_ELIGIBLE (neuro script only) |
| `params_motion_*.json` | Full parameter provenance record |
| `thresholds.tsv` | QC threshold documentation |
| `column_dictionary_motion.csv` | Column definitions for all output CSVs |

### Movement Group Comparison Outputs

| File | Description |
|---|---|
| `group_compare_summary.csv` | Statistical comparison (effect sizes, CIs, p-values) |
| `per_file_metrics.csv` | Per-file metric values for both groups |
| `*_box_swarm.png/pdf` | Box-swarm plots for each metric |
| `effects_forest_*.png/pdf` | Forest plot of Hedges' g effect sizes |
| `conclusions_panel_*.png/pdf` | Full statistical conclusions panel |

### Calcium Pipeline Outputs

**Step 0** — Motion Tracking outputs are listed above (shift CSVs, QC flags, `roi_eligible_list.txt`).

**Step 1** — ROI Detection (per file):

| File | Description |
|---|---|
| `*_rigid.tif` | Motion-corrected BigTIFF |
| `*_rois.mat` | Full ROI data structure (masks, traces, metrics) |
| `*_roi_metrics_v45.csv` | Per-ROI metrics (SNR, event rate, area, shape) |
| `*_analysis_v45.png/pdf` | 8-panel diagnostic figure |
| `*_roi_map_v45.png/pdf` | Clean ROI map with numbered outlines |
| `*_roi_gallery_v45.png/pdf` | Zoomed individual ROI thumbnails |
| `*_traces_v45.png/pdf` | ΔF/F calcium traces (paginated if >16 ROIs) |
| `column_dictionary_v45.txt` | Documentation of all CSV columns |
| `pipeline_summary_v45.csv` | Summary across all processed files |

**Step 2** — Quality Assessment:

| File | Description |
|---|---|
| `sample_quality_assessment.csv` | All samples with quality categories and metrics |
| `eligible_samples_*.txt` | Lists of usable samples per condition |
| `quality_comparison_summary.csv` | Summary statistics across conditions |
| `quality_diagnostic_figures/` | Per-sample diagnostic plots |

**Step 3** — ROI Refinement (per file):

| File | Description |
|---|---|
| `*_refined.mat` | Refined ROI data (same structure, fewer ROIs) |
| `*_refinement_log.csv` | Each removed ROI with the specific rejection reason |
| `*_rejection_examples.png/pdf` | Visualization showing why ROIs were removed |
| `refinement_summary.csv` | Cross-condition refinement statistics |

**Step 4** — Statistical Comparison:

| File | Description |
|---|---|
| `statistical_comparison.csv` | Full statistical results (effect sizes, CIs, p-values) |
| `forest_plot.png/pdf/svg` | Effect sizes with bootstrap confidence intervals |
| `removal_breakdown.png/pdf/svg` | ROI removal reasons by condition |
| `analysis_report.txt` | Comprehensive text report |

## Reproducing the Manuscript Results

To reproduce the results reported in the paper:

1. **Install dependencies** — MATLAB R2024b, NoRMCorre, Bio-Formats, and the three MATLAB toolboxes listed above.
2. **Set your imaging parameters** — Pixel size (266/1024 µm for our system) and frame rates (1 Hz for GFP, 33.33 Hz for calcium).
3. **Movement analysis (Workflow A):**
   1. Run `analyze_motion_and_QC.m` on each condition's GFP video folder
   2. Run `compare_motion_groups.m` for full-duration comparison
   3. Run `compare_motion_groups.m` again with `cfg.time_window_min = [0 30]` for the first-30-minute exploratory analysis
4. **Calcium imaging (Workflow B)** — run in order:
   1. `analyze_motion_and_QC_neuro.m` — produces shift CSVs and `roi_eligible_list.txt`
   2. `COMPLETE_PIPELINE_RUN_ME_v4_5.m` — reads shift CSVs, produces ROIs and traces
   3. `Assess_ROI_Quality.m` — flags artifacts and neuropil contamination
   4. `ROI_REFINEMENT_v3_2.m` — removes false-positive ROIs
   5. `COMPARE_NEURO_CONDITIONS_v3.m` — statistical comparison on refined data

All default parameters match those used in the manuscript.

## License

MIT License

Copyright (c) 2026 David Reynolds

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Citation

If you use this code, please cite:

> Reynolds, D., Artenyan, E., Nazaryan, H., Shanakian, E., Chen, E., Abramian, V., Ghashghaei, A., Sahabi, K., Safieh, F., Cazaubon, A., Momjian, N., Sunthorncharoenwong, J., Naddour, G., Arisaka, K. Multimodal immobilization of second-instar *Drosophila melanogaster* larvae using PF-127 hydrogel and diethyl ether for calcium imaging. *Journal of Neuroscience Methods* (submitted).

## Contact

For questions about the code or to request the original lab-internal versions, contact David Reynolds (davidrey555@gmail.com).
