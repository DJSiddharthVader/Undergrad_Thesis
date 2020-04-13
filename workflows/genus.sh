#!/bin/bash

#Globals
DATA_DIR="/home/sid/thesis_SidReed/data/"
SCRIPT_DIR="/home/sid/thesis_SidReed/scripts"

#Input Files
GENOMIC_DIR="${DATA_DIR}/genomic_data"
ALL_DATA="${DATA_DIR}/all_data.tsv"
NUC_DIR="./nucleotide"
PROT_DIR="./protein"
GENOME_DIR="./genome"

#Output Files/Dirs
MMSEQ_DIR="./mmseq_clustering"
FASTA_HEADER_DATA="./fasta_headers_info.json"
GENE_FAMILIES="./gene_families.json"
ALL_NUC_FASTA="./all_nucleotides.fna"
ALL_PROT_FASTA="./all_proteins.faa"
PA_MATRIX="./binary_pa_matrix.csv"
SPECIES_TREE_FILE="./species_tree_WGS/WGS_species_tree.con.tre"


function setUpFiles() {
    genus="$1"
    processes="$2"
    logfile="$3"
    mkdir -p "$genus" && cd "$genus"
    if [[ -f "$ALL_NUC_FASTA" ]]; then
        :
    else
        python "$SCRIPT_DIR"/get_genus_fastas.py "$ALL_DATA" "$genus" "$GENOMIC_DIR"
        python "$SCRIPT_DIR"/shortenFastaheaders.py "$NUC_DIR"
        cat "$NUC_DIR"/*.fna >| "$ALL_NUC_FASTA"
    fi
    if [[ -f "$ALL_PROT_FASTA" ]]; then
        :
    else
        python "$SCRIPT_DIR"/shortenFastaheaders.py "$PROT_DIR"
        cat "$PROT_DIR"/*.faa >| "$ALL_PROT_FASTA"
    fi
}
function clusterGenes() {
    genus="$1"
    processes="$2"
    logfile="$3"
    "$SCRIPT_DIR"/mmseqWrapper.sh "$ALL_PROT_FASTA" "$genus" "$MMSEQ_DIR"
    python "$SCRIPT_DIR"/parse_mmseq_clusters.py ./mmseq_clustering/clusters_"$genus".tsv "$ALL_PROT_FASTA"
    python "$SCRIPT_DIR"/createPAMatrix.py "$GENE_FAMILIES" "$PROT_DIR"

}
function buildTrees() {
    genus="$1"
    processes="$2"
    logfile="$3"
    python "$SCRIPT_DIR"/build_16S_species_tree_from_families.py "$ALL_NUC_FASTA" "$FASTA_HEADER_DATA"
    python "$SCRIPT_DIR"/build_WGS_species_tree.py "$genus" "$GENOME_DIR" $processes
    python "$SCRIPT_DIR"/generate_gene_trees.py "$GENE_FAMILIES" "$ALL_NUC_FASTA" $processes
}
function generateResults() {
    genus="$1"
    processes="$2"
    logfile="$3"
    Rscript "$SCRIPT_DIR"/network_builder.R ./ "$SPECIES_TREE_FILE" 1000 50 $processes
    Rscript "$SCRIPT_DIR"/markophylo_indel_estiamtes.R "$SPECIES_TREE_FILE" "$PA_MATRIX" "$ALL_DATA" 2>| markophylo.errs
    python "$SCRIPT_DIR"/parse_markophylo_output.py ./
}
function main() {
    genus="$1"
    processes="$2"
    logfile="$genus".log
    echo $genus
    setUpFiles $genus $processes $logfile
    clusterGenes $genus $processes $logfile
    buildTrees $genus $processes $logfile
    generateResults $genus $processes $logfile
}


main "$1" "$2"
#arg 1 is the genus, arg 2 is number of processes for the scripts that take that argument
#creates a directory in the current dir named $1 and puts all files results in that
