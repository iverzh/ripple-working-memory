#Ripple Analysis for Working Memory Task

This repository contains MATLAB code for analyzing ripple oscillations in human intracranial EEG data during a Sternberg working memory task. The analysis includes ripple detection, characterization, co-occurrence analysis across brain regions, and single-trial visualization.

## Data Availability

All data analyzed in this repository are from a previously published  **publicly available** through the DANDI Archive:

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

The dataset includes recordings from patients with pharmacologically intractable epilepsy who were implanted with depth electrodes for clinical monitoring. All data have been de-identified in accordance with HIPAA regulations.


## System Requirements

### Software Dependencies
- **MATLAB** R2019b or later (tested on R2019b)
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Image Processing Toolbox (for visualization)
  
### Required MATLAB Packages
The following packages should be available in your MATLAB path:
- **EEGLAB** (version 2022.0 or later) - for EEG data processing
- **boundedline** - for plotting mean ± SEM
- **violinplot** - for distribution visualization
- **slanCM** - for custom colormaps
- **brewermap** - for ColorBrewer colormaps

### Operating Systems
The code has been tested on:
- Linux (CentOS 7, Ubuntu 20.04)

### Hardware Requirements
- **Minimum**: 16 GB RAM, quad-core processor
- **Recommended**: 32+ GB RAM, 8+ core processor for processing large datasets
- **Storage**: Minimum 50 GB free space for data processing and output
- No specialized hardware (e.g., GPU) is required


## Installation Guide

### Step 1: Clone the Repository
```bash
git clone https://github.com/iverzh/ripple-working-memory.git
cd ripple-working-memory
```

### Step 2: Install Dependencies

#### EEGLAB
1. Download EEGLAB from [https://sccn.ucsd.edu/eeglab/download.php](https://sccn.ucsd.edu/eeglab/download.php)
2. Extract to a directory (e.g., `/path/to/eeglab2022.0/`)
3. Add to MATLAB path:
   ```matlab
   addpath(genpath('/path/to/eeglab2022.0/'))
   ```

#### matnwb (for reading DANDI NWB files)
```matlab
% Clone and setup matnwb
% From command line:
git clone https://github.com/NeurodataWithoutBorders/matnwb.git
cd matnwb
git checkout v2.6.0  % or latest stable version

% In MATLAB:
cd('path/to/matnwb')
addpath(genpath(pwd));
generateCore();  % Generate MATLAB classes for NWB
```

#### Plotting Functions
Install the required plotting functions from MATLAB File Exchange or GitHub:

```matlab
% Install boundedline
% https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

% Install violinplot
% https://github.com/bastibe/Violinplot-Matlab

% Install brewermap
% https://github.com/DrosteEffect/BrewerMap

% Install slanCM
% https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap
```

### Step 3: Update Path Variables
Edit the path variables at the beginning of each script to match your local directory structure:

```matlab
% In ripple_characteristics_cleaned.m
dataDirectoryOrig = '/path/to/your/data/preprocess/OrigUpload';
dataDirectory = '/path/to/your/data/preprocess/data_1kHz';
matExportFolder = '/path/to/your/output/matFiles';
exportDir = '/path/to/your/output/processedResults';
exportDirFigs = '/path/to/your/output/figures';
```

### Installation Time
Typical installation time: **1 minute** on a standard desktop computer with a stable internet connection.
Entire working memory dataset is **22 GB**

## Demo

### Downloading Sample Data from DANDI

Before running the analysis, you'll need to download data from the DANDI Archive:

#### Using the DANDI CLI

1. **Install the DANDI CLI**:
   ```bash
   pip install dandi
   ```

2. **Download the dataset**:
   ```bash
   # Download the entire dataset (Warning: this may be large!)
   dandi download https://dandiarchive.org/dandiset/000673
   
   # Or download specific subjects only
   dandi download https://dandiarchive.org/dandiset/000673 -i "sub-01/*"
   ```

3. **Access the data**: The data will be downloaded in NWB format. You can read NWB files in MATLAB using [matnwb](https://github.com/NeurodataWithoutBorders/matnwb):

### Data Preprocessing

The analysis scripts expect preprocessed data in MATLAB `.mat` format. The DANDI dataset provides raw NWB files, which need to be converted:

1. **Extract LFP data** from NWB files and save as MATLAB matrices
2. **Extract unit data** (spike times, waveforms, classifications)
3. **Extract task events** (stimulus timing, behavioral responses)
4. **Run ripple detection** using the detection pipeline (see [ripple-detection](https://github.com/iverzh/ripple-detection) repository)



### Expected Run Time
- **Demo run time**: 5-10 minutes per subject on a standard desktop (Intel i7, 16GB RAM)
- **Full dataset**: 2-4 hours for 20-30 subjects
- The bottleneck is primarily time-frequency decomposition (ERSP calculation)

## Instructions for Use

### Data Format Requirements

The code expects data in the following format. **Note**: The source data from DANDI (dandiset 000673) are in NWB format and need to be converted to these MATLAB formats before analysis.

#### Converting from DANDI/NWB Format

The DANDI dataset (000673) provides all necessary data in NWB format:
- **LFP recordings**: Available in `nwb.acquisition`
- **Single-unit data**: Available in `nwb.units` (includes spike times, waveforms, quality metrics)
- **Electrode locations**: Available in `nwb.electrodes`
- **Task events**: Available in `nwb.intervals` or `nwb.trials`

You'll need to extract these from NWB files and convert to the MATLAB formats described below. See the [matnwb documentation](https://neurodatawithoutborders.github.io/matnwb/) for detailed extraction instructions.

#### 1. LFP Data Files
```matlab
% File: <subject>_LFP_micro.mat
% Variables:
%   - data: [channels × time] matrix of LFP voltage (?V)
%   - times: [1 × time] vector of timestamps (seconds)
%   - chan_locations: cell array of channel location labels
%   - fs: sampling frequency (Hz)
```

#### 2. Ripple Statistics Files
Output from [ripple detection pipeline](https://github.com/iverzh/ripple-detection)
```matlab
% File: <subject>_ripple_stats_wake_NC_1kHz_template_z25.mat
% Variables in rippleStats structure:
%   - locs: {1 × channels} cell array of ripple peak indices
%   - window: {1 × channels} cell array of [start, end] times
%   - duration: {1 × channels} cell array of ripple durations (ms)
%   - oscFreq: {1 × channels} cell array of oscillation frequencies (Hz)
%   - rippleAmp: {1 × channels} cell array of ripple amplitudes (uV)
%   - density: {1 × channels} cell array of ripple density (events/min)
%   - chanLabels: {1 × channels} cell array of channel labels
%   - recordingLength: total recording length (samples)
%   - fs: sampling frequency (Hz)
```

#### 3. Unit (Spike) Data Files
```matlab
% File: <subject>_unit.mat
% Format: Cell array with each row containing:
%   Column 1: Channel number
%   Column 2: Spike times (seconds)
%   Column 3: Unit type ('pyr', 'int', 'mult')
%   Column 4: Cluster ID
%   Column 5: Brain region label
%   Column 6: Trough-to-peak width (samples)
%   Column 7: Burst index
```

#### 4. Task Data Files
```matlab
% File: <subject>_task.mat
% Variables:
%   - start_time: [trials × 1] trial start times (seconds)
%   - loads: [trials × 1] memory load for each trial (1-3)
%   - timestamps_Encoding1/2/3: Image presentation times
%   - timestamps_Maintenance: Maintenance period start
%   - timestamps_Probe: Probe presentation time
%   - timestamps_Response: Response time
%   - PicIDs_Encoding1/2/3: Image IDs for encoding
%   - PicIDs_Probe: Probe image ID
%   - probe_in_out: [trials × 1] logical, true if probe was in encoding set
```

### Running the Analysis

#### 1. Main Ripple Characteristics Analysis
```matlab
% This script analyzes ripple properties across brain regions
ripple_characteristics_cleaned

% Key sections:
% - Lines 1-100: Setup and parameter definition
% - Lines 101-500: Subject processing loop
% - Lines 501-end: Visualization and statistics
```

**Customize analysis parameters:**
```matlab
% In the script, modify these parameters:
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};  % Regions to analyze
win = 1000;  % CCG window size (ms)
binSize = 1;  % Histogram bin size (ms)
```

#### 2. Single-Trial Co-Ripple Visualization
```matlab
% This script generates single-trial plots of co-ripple events
ripple_single_trial_cleaned

% Sections:
% - Section 1: Co-ripple events during successful retrieval
% - Section 2: Single ripple events during unsuccessful retrieval
```

**Customize visualization:**
```matlab
% Modify these parameters:
ripWin = 200;  % Window size for visualization (ms)
fs = 1e3;  % Sampling frequency (Hz)
```

### Output Files

The analysis generates several types of output:

1. **Figures** (saved as PDF in `exportDirFigs/`):
   - Ripple waveforms aligned to ripple center
   - Distribution of ripple metrics (density, duration, amplitude, frequency)
   - Channel and unit count distributions
   - Time-frequency representations (ERSP)
   - Co-ripple rate analyses

2. **Statistics** (printed to console):
   - Mean ± SEM for ripple metrics by region
   - Channel and unit counts by region
   - Co-ripple occurrence rates

3. **Single-trial visualizations** (if using `ripple_single_trial_cleaned.m`):
   - Individual trial plots showing LFP, filtered ripple band, and spike times
   - Comparison of encoding and retrieval periods

### Processing Your Own Data

To adapt this code for your own dataset:

1. **Format your data** according to the specifications above
2. **Modify file paths** in the script headers
3. **Update region labels** to match your anatomical parcellation:
   ```matlab
   % Update these mapping functions:
   locations = strrep(locations, 'your_region_name', 'STANDARD_NAME');
   ```
4. **Adjust filtering parameters** if your ripple band differs:
   ```matlab
   % Modify bandpass filter
   [b,a] = butter(3, [70 100]/(fs/2));  % Default: 70-100 Hz
   ```
5. **Run the analysis** and check outputs for quality


## Citation

If you use this code in your research, please cite:

**For this analysis code:**
```
Verzhbinsky, I. A., Daume, J., Rutishauser, U., & Halgren, E. (2025). 
Cross-region neuron co-firing mediated by ripple oscillations supports 
distributed working memory representations. bioRxiv, 2025-09.
```

**For the source dataset (DANDI 000673):**
```
Daume, J., Kami?ski, J., Schjetnan, A. G. P., Salimpour, Y., Khan, U., 
Kalia, S. K., Valiante, T. A., Anderson, A., Mamelak, A. N., & Rutishauser, U. 
Control of working memory by phase-amplitude coupling of human hippocampal neurons. 
Nature (2024). https://doi.org/10.1038/s41586-024-07309-z

Dataset available at: https://dandiarchive.org/dandiset/000673
```


**For the ripple detection algorithm:**
```
[Ripple detection paper citation - if using iverzh/ripple-detection]
```

## Related Resources

- **Source Dataset**: https://dandiarchive.org/dandiset/000673
- **Ripple Detection Code**: https://github.com/iverzh/ripple-detection
- **DANDI Archive**: https://dandiarchive.org
- **DANDI CLI**: https://github.com/dandi/dandi-cli
- **matnwb (MATLAB NWB interface)**: https://github.com/NeurodataWithoutBorders/matnwb
- **NWB Documentation**: https://nwb-overview.readthedocs.io/

