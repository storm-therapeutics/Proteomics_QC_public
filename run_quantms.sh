#!/bin/bash

FASTA=/home/hendrik.weisser/Data/Proteomics/SST/Human_reference_proteome_TD.fasta
# experimental design file for mzML/raw input:
EXP_DESIGN_MZML=/home/hendrik.weisser/Data/Proteomics/SST/QuantMS/single_runs/exp_design_mzML.tsv
EXP_DESIGN_RAW=/home/hendrik.weisser/Data/Proteomics/SST/QuantMS/single_runs/exp_design_raw.tsv

if [ $# -eq 0 ]; then
    echo "Usage: run_quantms.sh OUT_DIR FILE..."
fi

outdir=$1
shift 1 # argument list now starts at the 2nd argument

if [ $# -eq 0 ]; then
    echo "Error: no input files specified"
    exit 1
fi

# get abs. paths to input files before switching directory:
paths=`realpath $@`

mkdir -p $outdir
cd $outdir

for p in $paths; do
    f=`basename $p`
    echo "Working on" $f
    ext="${f##*.}" # file extension
    ext_lower=${ext,,} # convert to lowercase
    if [ $ext_lower = "mzml" ]; then
        EXP_DESIGN=$EXP_DESIGN_MZML
        ln -f -s $p current.mzML # overwrite if destination exists
    elif [ $ext_lower = "raw" ]; then
        EXP_DESIGN=$EXP_DESIGN_RAW
        ln -f -s $p current.raw # overwrite if destination exists
    else
        echo "Error: unsupported file type:" $ext
        exit 1
    fi
    nextflow run nf-core/quantms -r dev -profile docker --input $EXP_DESIGN --database $FASTA --search_engines comet --variable_mods "Oxidation (M),Acetyl (Protein N-term)" --max_precursor_charge 5 --min_peptide_length 7 --FDR_level psm-level-fdrs --psm_pep_fdr_cutoff 0.01 --protein_level_fdr_cutoff 0.01 --max_memory 48.GB --skip_post_msstats --outdir results_${f%.*} --labelling_type "label free sample" --acquisition_method dda
done

# clean up symlinks:
rm -f current.mzML current.raw
