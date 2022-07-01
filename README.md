# Proteomics_QC
Code for proteomics quality control -- processing LC-MS/MS data from "system suitability test" (SST) samples

## General idea
The same kind of sample (HeLa lysate) is analysed by LC-MS/MS before every relevant experiment to monitor the function of the mass spectrometer (and liquid chromatography system).
If the resulting data is processed in a consistent way, instrument performance can be tracked over time.
Data analysis should be quick and convenient to use, and generate a report with helpful statistics and visualisations.

## Data analysis workflow
SST samples are processed one-at-a-time using the ["nf-core/quantms"](https://nf-co.re/quantms/dev) Nextflow workflow using suitable settings.
The main result is an [mzTab](https://github.com/HUPO-PSI/mzTab) file containing identification and quantification information on PSM, peptide and protein level.
mzTab files from multiple SST runs can be flexibly merged and a QC report in HTML format generated from the combined results.

## Scripts

### 1. `run_quantms.sh`
Shell script for running the "nf-core/quantms" workflow on multiple input files (.raw or .mzML), but processing each file separately.  
Usage: `./run_quantms.sh out_dir /path/to/in1.raw /path/to/in2.raw [...]`  
This will create the `out_dir` directory and process each input file, storing results in subdirectories named after the input files (`out_dir/results_in1` etc.).

### 2. `proteomics_qc.R`
Executable R script for merging mzTab files from "nf-core/quantms" results and generating a QC report (HTML).  
Usage:  
```
1. Rscript proteomics_qc.R results_in1 results_in2 [...]
2. Rscript proteomics_qc.R merged.mzTab results_in1 results_in2 [...]
3. Rscript proteomics_qc.R merged.mzTab
```
The first call will read results from all input directories and create a merged mzTab files (`merged.mzTab`) and QC report (`qc_report.html`) in the current directory.  
The second call does the same, but merges the results with an existing mzTab file from a previous invocation of the script.
This saves processing time when new samples are compared to older ones.  
The third call recreates the QC report from an existing mzTab file.
This is useful when the report format has changed.

### 3. `qc_report.Rmd`
[R Markdown](https://bookdown.org/yihui/rmarkdown/) file containing code and text for the QC report.
Using a (merged) mzTab file for input data, this gets rendered to an HTML file by the `proteomics_qc.R` script.
