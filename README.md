# Ripple Analysis for Working Memory Task

This repository contains MATLAB and Python code accompanying:

> Verzhbinsky, I. A., Daume, J., Rutishauser, U., & Halgren, E.
> *Cross-region neuron co-firing mediated by ripple oscillations supports
> distributed working memory representations.* Nature Neuroscience (in press, 2026).

The code analyzes ripple oscillations and single-unit activity recorded from
human intracranial electrodes during a Sternberg working memory task. It
covers preprocessing of the public DANDI dataset, ripple characterization,
cross-regional co-rippling and co-firing, ripple-locked replay of encoding
co-firing patterns, and single-trial / anatomical visualization.

## Data Availability

All data analyzed in this repository are **publicly available** through the
DANDI Archive:

**DANDI Dataset 000673**: *Control of working memory by phase amplitude coupling of human hippocampal neurons*
- **URL**: https://dandiarchive.org/dandiset/000673
- **Contents**:
  - Local field potential (LFP) recordings from human intracranial electrodes
  - Single-unit spike times and waveforms
  - Electrode coordinates and anatomical locations
  - Behavioral task data (Sternberg working memory paradigm)
  - Experimental parameters and anonymized patient metadata
- **Format**: Neurodata Without Borders (NWB) format
- **License**: Creative Commons (see DANDI dataset for details)

The dataset includes recordings from patients with pharmacologically
intractable epilepsy who were implanted with depth electrodes for clinical
monitoring. All data have been de-identified in accordance with HIPAA
regulations.

## Repository Structure

```
ripple-working-memory/
├── LICENSE
├── README.md
└── code/
    ├── preprocess_LFP_units.ipynb        # NWB → .mat preprocessing (Python)
    ├── RippleTaskAnalysis_WM.m           # Ripple/firing rates across task phases
    ├── coripple_cofire_analysis.m        # Cross-region co-firing during co-ripples
    ├── coripple_unit_replay.m            # Encoding → probe replay analysis
    ├── coripple_unit_replay_timecourse.m # PSTH of replay around ripple peaks
    ├── plot_ripples.m                    # Ripple waveforms, ERSP, distributions
    ├── plotReplaySingleTrials.m          # Single-trial co-ripple/replay visualizations
    ├── plotMNI.m                         # Electrode locations on MNI brain surface
    └── util/                             # Helpers + bundled plotting libs
        ├── aggregateRegionData.m
        ├── aggregateRegionDataShuff.m
        ├── aggregateRegionDataSubj.m
        ├── compute_sttc.m                # Spike-time tiling coefficient
        ├── LoadSpikeTimes.m
        ├── chiSquaredTest.m, cohens_d.m, fdr_bh.m,
        ├── pairedPermutationTest.m, permutationTestMed.m,
        ├── computeSEM.m, bounds2mask.m, mask2bounds.m,
        ├── exportToVTK.m, savepdf.m,
        ├── bar_striped.m, errobarPatch.m, hline.m, vline.m, tight_subplot.m,
        ├── brewermap.m
        └── boundedline/  dabarplot/  daboxplot/  daviolinplot/  slanCM/  subtightplot/
```

## System Requirements

### Software Dependencies
- **MATLAB** R2019b or later (tested on R2019b and R2023a)
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Image Processing Toolbox (for visualization)
- **Python** 3.9+ with `pynwb`, `numpy`, `scipy`, `h5py`, and `jupyter`
  (only required for `preprocess_LFP_units.ipynb`)

### Bundled MATLAB Packages
The following plotting/utility packages are **bundled in `code/util/`** and
do not need to be installed separately — just add `code/util/` (recursively)
to your MATLAB path:
- **boundedline** — mean ± SEM line plots
- **dabarplot / daboxplot / daviolinplot** — distribution visualization
- **slanCM** and **brewermap** — colormaps
- **subtightplot / tight_subplot** — multi-panel layouts
- **fdr_bh** — Benjamini–Hochberg FDR correction

### Optional / External
- **EEGLAB** (2022.0 or later) — only needed if you re-run the upstream
  ripple-detection pipeline; not required for the analyses in this repo.
- **matnwb** — only needed if you prefer reading NWB files in MATLAB rather
  than using the supplied Python preprocessing notebook.

### Operating Systems
The code has been tested on:
- Linux (CentOS 7, Ubuntu 20.04 / 22.04)
- macOS 13–14

### Hardware Requirements
- **Minimum**: 16 GB RAM, quad-core processor
- **Recommended**: 32+ GB RAM, 8+ core processor for processing large datasets
- **Storage**: ~50 GB free space (full DANDI dataset is ~22 GB)
- No specialized hardware (e.g., GPU) is required

## Installation Guide

### Step 1: Clone the Repository
```bash
git clone https://github.com/iverzh/ripple-working-memory.git
cd ripple-working-memory
```

### Step 2: Set Up MATLAB Path
```matlab
addpath('code');
addpath(genpath('code/util'));   % adds bundled plotting libraries
```

### Step 3: Set Up Python Environment (for preprocessing only)
```bash
pip install pynwb numpy scipy h5py jupyter
```

### Step 4: Update Path Variables
Each analysis script defines its data and output paths in a header block
near the top of the file. Edit these to match your local layout, e.g.:

```matlab
dataDirectory   = '/path/to/preprocessed/data_1kHz';
matExportFolder = '/path/to/output/matFiles';
exportDir       = '/path/to/output/processedResults';
exportDirFigs   = '/path/to/output/figures';
```

### Installation Time
Typical installation time: **~1 minute** on a standard desktop with a
stable internet connection (excluding DANDI dataset download).

## Demo

### Step 1 — Download the DANDI Dataset
```bash
pip install dandi
dandi download https://dandiarchive.org/dandiset/000673
```
You can restrict to a single subject for a quick demo:
```bash
dandi download https://dandiarchive.org/dandiset/000673 -i "sub-01/*"
```

### Step 2 — Convert NWB → MATLAB
Open and run `code/preprocess_LFP_units.ipynb`. Update the `dandi_root`
and `out_root` paths in the first cell. For each session it produces:

- `<subject>_LFP_micro.mat` — LFP traces, channel labels, MNI coordinates,
  trial epoch boundaries, sampling rate, and timestamps
- `<subject>_unit.mat` — single-unit spike times and 1-indexed electrode IDs
- `<subject>_task.mat` — Sternberg trial events, loads, and probe accuracy

The notebook applies a 10-second edge pad around each trial window to
minimize filter edge artifacts during downstream ripple detection.

### Step 3 — Detect Ripples
Ripple detection is implemented in a separate repository:
[iverzh/ripple-detection](https://github.com/iverzh/ripple-detection).
Run that pipeline on the `_LFP_micro.mat` files to produce:

- `<subject>_ripple_stats_wake_NC_1kHz_template_z25.mat`

### Step 4 — Run the Analyses
With paths configured, scripts can be run independently in MATLAB. A
typical order is:

```matlab
RippleTaskAnalysis_WM            % task-evoked ripple/firing rates per region
plot_ripples                     % mean ripple waveforms, ERSPs, distributions
coripple_cofire_analysis         % cross-region co-firing during co-ripples
coripple_unit_replay             % encoding-pattern replay at probe
coripple_unit_replay_timecourse  % replay PSTH around ripple peaks
plotReplaySingleTrials           % single-trial replay/co-ripple figures
plotMNI                          % electrode coverage on MNI surface
```

### Expected Run Time
- **Per-subject demo**: 5–15 minutes on a standard desktop (Intel i7, 16 GB RAM)
- **Full dataset (≈30 subjects)**: 2–6 hours depending on which analyses are run
- The bottleneck is usually time-frequency decomposition (ERSP) and the
  pairwise co-firing computation in `coripple_cofire_analysis.m`.

## Instructions for Use

### Data Format Requirements

The MATLAB analysis scripts expect the following files per subject (all
produced by `preprocess_LFP_units.ipynb` plus the ripple-detection
pipeline).

#### 1. LFP Data Files
```matlab
% File: <subject>_LFP_micro.mat
% Variables:
%   - data:           [channels × time] matrix of LFP voltage (µV)
%   - times:          [1 × time] vector of timestamps (seconds)
%   - chan_locations: cell array of anatomical channel labels
%   - chan_xyz:       [channels × 3] MNI coordinates
%   - fs:             sampling frequency (Hz)
%   - trial_epochs:   [trials × 2] padded trial windows (samples)
```

#### 2. Ripple Statistics Files
Output of [ripple-detection](https://github.com/iverzh/ripple-detection):
```matlab
% File: <subject>_ripple_stats_wake_NC_1kHz_template_z25.mat
% rippleStats struct:
%   - locs:            {1 × channels} ripple peak indices
%   - window:          {1 × channels} [start, end] sample indices
%   - duration:        {1 × channels} ripple durations (ms)
%   - oscFreq:         {1 × channels} oscillation frequencies (Hz)
%   - rippleAmp:       {1 × channels} ripple amplitudes (µV)
%   - density:         {1 × channels} ripple density (events/min)
%   - chanLabels:      {1 × channels} channel labels
%   - recordingLength: total recording length (samples)
%   - fs:              sampling frequency (Hz)
```

#### 3. Unit (Spike) Data Files
```matlab
% File: <subject>_unit.mat
% Cell array, one row per unit:
%   Col 1: Channel number (1-indexed for MATLAB)
%   Col 2: Spike times (seconds)
%   Col 3: Unit type ('pyr', 'int', 'mult')
%   Col 4: Cluster ID
%   Col 5: Brain region label
%   Col 6: Trough-to-peak width (samples)
%   Col 7: Burst index
```

#### 4. Task Data Files
```matlab
% File: <subject>_task.mat
%   - start_time:               [trials × 1] trial start times (s)
%   - loads:                    [trials × 1] memory load (1 or 3)
%   - timestamps_Encoding1/2/3: image presentation times
%   - timestamps_Maintenance:   maintenance period start
%   - timestamps_Probe:         probe presentation time
%   - timestamps_Response:      response time
%   - PicIDs_Encoding1/2/3:     image IDs for encoding
%   - PicIDs_Probe:             probe image ID
%   - probe_in_out:             [trials × 1] true if probe was in encoding set
```

### Script Reference

| Script | Purpose |
| --- | --- |
| `preprocess_LFP_units.ipynb` | Convert DANDI NWB sessions into the `_LFP_micro.mat`, `_unit.mat`, and `_task.mat` files used by all MATLAB scripts. |
| `RippleTaskAnalysis_WM.m` | Ripple rate, firing rate, and ripple–spike coupling across baseline / encoding / maintenance / probe phases and load conditions; LME and FDR-corrected statistics across regions (OFC, ACC, SMA, AMY, HIP). |
| `plot_ripples.m` | Mean ripple waveforms (with SEM), ERSP, violin plots of ripple density / duration / amplitude / frequency, channel and unit count summaries, co-ripple rate vs. tract length. |
| `coripple_cofire_analysis.m` | Pairwise spike co-firing in a 25 ms window across regions, conditioned on co-ripple vs. non-ripple periods, separated by task phase and load. Includes anatomical-distance correlation and load-dependent modulation. |
| `coripple_unit_replay.m` | Tests whether co-firing patterns at encoding are reinstated at probe, and whether reinstatement is enriched during co-ripples. Outputs subject-level replay rates and RT correlations. |
| `coripple_unit_replay_timecourse.m` | PSTH of replay events around ripple peaks, separating co-ripple, single-region ripple, and shuffled controls; computed for both encoding and probe periods. |
| `plotReplaySingleTrials.m` | Four-panel single-trial figures: raw LFP, 70–100 Hz band, spike rasters, and ripple windows for two regions during encoding and probe. |
| `plotMNI.m` | Builds an MNI cortical surface from a T2 template, overlays FreeSurfer pial surfaces, and exports per-region electrode point clouds as VTK. |

### Customizing Analysis Parameters

Common knobs near the top of each script:

```matlab
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};   % regions to include
fs      = 1e3;                                    % sampling frequency (Hz)
win     = 1000;                                   % CCG / PSTH half-window (ms)
binSize = 1;                                      % histogram bin size (ms)
ripWin  = 200;                                    % single-trial plot window (ms)
[b,a]   = butter(3, [70 100] / (fs/2));           % ripple bandpass
```

### Output Files

The analysis scripts produce:

1. **Figures** (PDFs in `exportDirFigs/`): ripple waveforms, ERSPs,
   distribution violins, electrode/unit count pies, ripple-rate timecourses,
   co-ripple and co-firing heatmaps, replay rate bar/box plots, single-trial
   panels, MNI electrode renders.
2. **Statistics** (printed to console + saved `.mat`): means ± SEM by region,
   linear mixed-effects model coefficients, Wilcoxon and permutation
   p-values with FDR correction, chi-square tests for replay enrichment.
3. **Intermediate aggregates** (in `exportDir/`): per-subject and per-region
   `.mat` files written by `aggregateRegionData*` so plotting scripts can
   re-run without recomputing the heavy analyses.

## Citation

**For this analysis code and the accompanying paper:**
```
Verzhbinsky, I. A., Daume, J., Rutishauser, U., & Halgren, E. (2026).
Cross-region neuron co-firing mediated by ripple oscillations supports
distributed working memory representations. Nature Neuroscience (in press).

Preprint: https://www.biorxiv.org/content/10.1101/2025.09.04.674061v1
```

**For the source dataset (DANDI 000673):**
```
Daume, J., Kamiński, J., Schjetnan, A. G. P., Salimpour, Y., Khan, U.,
Kalia, S. K., Valiante, T. A., Anderson, A., Mamelak, A. N., & Rutishauser, U.
Control of working memory by phase-amplitude coupling of human hippocampal neurons.
Nature 629, 393–401 (2024). https://doi.org/10.1038/s41586-024-07309-z

Dataset: https://dandiarchive.org/dandiset/000673
```

**For the ripple detection algorithm:**
```
Dickey, C. W., Verzhbinsky, I. A., Jiang, X., Rosen, B. Q., Kajfez, S.,
Stedelin, B., Shih, J. J., Ben-Haim, S., Raslan, A. M., Eskandar, E. N.,
Gonzalez-Martinez, J., Cash, S. S., & Halgren, E.
Widespread ripples synchronize human cortical activity during sleep,
waking, and memory recall.
Proceedings of the National Academy of Sciences 119(28), e2107797119 (2022).
https://doi.org/10.1073/pnas.2107797119

Code: https://github.com/iverzh/ripple-detection
```

## Related Resources

- **Source Dataset**: https://dandiarchive.org/dandiset/000673
- **Ripple Detection Code**: https://github.com/iverzh/ripple-detection
- **DANDI Archive**: https://dandiarchive.org
- **DANDI CLI**: https://github.com/dandi/dandi-cli
- **pynwb (Python NWB interface)**: https://pynwb.readthedocs.io
- **matnwb (MATLAB NWB interface)**: https://github.com/NeurodataWithoutBorders/matnwb
- **NWB Documentation**: https://nwb-overview.readthedocs.io/

## License

This code is released under the MIT License (see `LICENSE`).
