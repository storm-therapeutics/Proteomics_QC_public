#!/bin/bash

FASTA=/home/hendrik.weisser/Data/Proteomics/SST/Human_reference_proteome_TD.fasta
# experimental design file (template):
EXP_DESIGN=/home/hendrik.weisser/Data/Proteomics/SST/QuantMS/exp_design_single.tsv

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
    # replace path to input file in exp. design template:
    sed s@PATH@"$p"@ $EXP_DESIGN > exp_design.tsv

    nextflow run nf-core/quantms -r dev -profile docker --input exp_design.tsv --database $FASTA --search_engines comet --variable_mods "Oxidation (M),Acetyl (Protein N-term)" --max_precursor_charge 5 --min_peptide_length 7 --FDR_level psm-level-fdrs --psm_pep_fdr_cutoff 0.01 --protein_level_fdr_cutoff 0.01 --max_memory 48.GB --skip_post_msstats --outdir results_${f%.*} --labelling_type "label free sample" --acquisition_method dda
done
