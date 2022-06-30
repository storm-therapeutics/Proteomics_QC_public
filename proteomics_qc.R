#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=FALSE)
## find path to input script (for "source" call later):
pattern <- "^--file="
script.path <- sub(pattern, "", grep(pattern, args, value=TRUE))
script.dir <- normalizePath(dirname(script.path))

index <- which(args == "--args")
if (length(index) == 0) # no trailing arguments
{
  cat("Usage: Rscript proteomics_qc.R [merged.mzTab] quantms_dir...\n",
      "Generate a proteomics QC report from 'nf-core/quantms' (Nextflow workflow) results.",
      "Inputs: List of 'quantms' result directories; optionally a merged mzTab file from a previous invocation.",
      "Outputs (created in the current directory): HTML report ('qc_report.html'), merged mzTab file ('merged.mzTab').", sep="\n")
  q("no")
}

args <- args[-(1:index)] # keep only trailing arguments
## for testing:
## args <- grep("results_", list.dirs("/home/hendrik.weisser/Data/Proteomics/SST/QuantMS/single_runs", recursive=FALSE), value=TRUE)

## optional merged mzTab file given?
if (grepl("\\.mztab$", args[1], ignore.case=TRUE)) {
  first.mztab <- args[1]
  if (!file.exists(first.mztab))
    stop("File not found: ", first.mztab)
  args <- args[-1]
} else {
  first.mztab <- NULL
}

for (arg in args) {
  if (!dir.exists(arg))
    stop("No such directory: ", arg)
}

if (length(args) > 0) {
  cat("Loading functions...\n")
  suppressMessages(source(file.path(script.dir, "merge_mzTab.R")))

  samples <- sub("^results_", "", basename(args))

  ## count number of MS2 spectra in each input file:
  ## (relevant line in mzML: <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" />)
  mzml.dir <- "mzmlindexing/out"
  n.ms2 <- sapply(args, function(dir) system2("grep", c("-c", "MS:1000580", file.path(dir, mzml.dir, "*.mzML")),
                                              stdout=TRUE))
  n.ms2 <- as.integer(n.ms2)

  cat("Merging mzTab data...\n")
  mztab.paths <- file.path(args, "proteomicslfq/out.mzTab")
  merged <- merge.mztab.files(mztab.paths)
  merged <- annotate.ms.runs(merged, n.ms2, paste0(samples, ".mzML"))

  if (!is.null(first.mztab)) {
    cat("Merging with previous mzTab file...\n")
    first <- load.mztab.file(first.mztab)
    merged <- merge.mztab.data(first, merged)
  }

  cat("Writing merged mzTab file...\n")
  writeMzTabData(merged, "merged.mzTab")
  merged.path <- normalizePath("./merged.mzTab")
} else { # just create the report from the given mzTab file
  merged.path <- normalizePath(first.mztab)
}

cat("Generating HTML report...\n")
library(rmarkdown)
render("/home/hendrik.weisser/Code/Proteomics_QC/qc_report.Rmd", output_dir=".",
       params=list(data=merged.path))
