#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=FALSE)
## find path to input script (for "source" call later):
pattern <- "^--file="
script.path <- sub(pattern, "", grep(pattern, args, value=TRUE))
script.dir <- normalizePath(dirname(script.path))
## keep only trailing arguments:
args <- args[-(1:which(args == "--args"))]

if (length(args) == 0)
  stop("List of directories (containing QuantMS results) required as arguments.")

for (arg in args) {
  if (!dir.exists(arg))
    stop("No such directory: ", arg)
}

## cat("Input directories:", args, sep="\n")

## for testing:
## args <- grep("results_", list.dirs("/home/hendrik.weisser/Data/Proteomics/SST/QuantMS/single_runs", recursive=FALSE), value=TRUE)
samples <- sub("^results_", "", basename(args))

## count number of MS2 spectra in each input file:
mzml.dir <- "mzmlindexing/out"
n.ms2 <- sapply(args, function(dir) system2("grep", c("-c", "MS:1000580", file.path(dir, mzml.dir, "*.mzML")),
                                            stdout=TRUE))
n.ms2 <- as.integer(n.ms2)

cat("Loading functions...\n")
suppressMessages(source(file.path(script.dir, "merge_mzTab.R")))

cat("Merging mzTab data...\n")
mztab.paths <- file.path(args, "proteomicslfq/out.mzTab")
merged <- merge.mztab.files(mztab.paths)
merged <- annotate.ms.runs(merged, n.ms2, paste0(samples, ".mzML"))

cat("Writing merged mzTab file...\n")
writeMzTabData(merged, "merged.mzTab")

cat("Generating HTML report...\n")
library(rmarkdown)
merged.path <- normalizePath("./merged.mzTab")
render("/home/hendrik.weisser/Code/Proteomics_QC/qc_report.Rmd", output_dir=".",
       params=list(data=merged.path))
