#!/usr/bin/env nextflow

//Globals
params.DATA_DIR="/home/sid/thesis_SidReed/data/"
params.SCRIPT_DIR="/home/sid/thesis_SidReed/scripts"
params.PROCESSES=8
params.logdir="/home/sid/thesis_SidReed/data/logs"

//Input Files
params.GENOMIC_DIR="${DATA_DIR}/genomic_data"
params.ALL_DATA="${DATA_DIR}/all_data.tsv"
params.NUC_DIR="nucleotide"
params.PROT_DIR="protein"
params.GENOME_DIR="genome"

//Output Files/Dirs
params.MMSEQ_DIR="mmseq_clustering"
params.FASTA_HEADER_DATA="fasta_headers_info.json"
params.GENE_FAMILIES="gene_families.json"
params.ALL_NUC_FASTA="all_nucleotides.fna"
params.ALL_PROT_FASTA="all_proteins.faa"
params.PA_MATRIX="binary_pa_matrix.csv"
params.SPECIES_TREE_FILE="species_tree_WGS/WGS_species_tree.con.tree"


Channel.fromPath(params.genera)
       .splitText()
       .map{file(it.trin())}
       .set{genera}


process main {
    input:
    params.SCRIPT_DIR
    params.PROCESSES
    params.DATA_DIR
    params.GENOMIC_DIR
    params.ALL_DATA
    params.NUC_DIR
    params.PROT_DIR
    params.GENOME_DIR
    params.MMSEQ_DIR
    params.FASTA_HEADER_DATA
    params.GENE_FAMILIES
    params.ALL_NUC_FASTA
    params.ALL_PROT_FASTA
    params.PA_MATRIX
    params.SPECIES_TREE_FILE
    val genus from genera

    output:
    file '$logdir/!{genus}.log' into finished

    """
    mkdir -p !{input} && cd !{input}
    python !{SCRIPT_DIR}/get_genus_fastas.py !{ALL_DATA} ./ !{input} !{GENOMIC_DIR}
    python !{SCRIPT_DIR}/shortenFastaheaders.py !{NUC_DIR}
    cat !{NUC_DIR}/*.fna >| !{ALL_NUC_FASTA}
    python !{SCRIPT_DIR}/shortenFastaheaders.py !{PROT_DIR}
    cat !{PROT_DIR}/*.faa >| !{ALL_PROT_FASTA}
    echo "single link clustering to gene families"
    !{scripts}/mmseqWrapper.sh all_proteins.faa !{genus} !{mmseq_dir}
    python !{SCRIPT_DIR}/parse_mmseq_clusters.py ./mmseq_clustering/clusters_!{input}.tsv !{ALL_PROT_FASTA}
    echo "creating P/A matrix"
    python !{SCRIPT_DIR}/createPAMatrix.py !{GENE_FAMILIES} !{PROT_DIR}
    echo "building species trees"
    // python !{SCRIPT_DIR}/build_16S_species_tree_from_families.py !{ALL_NUC_FASTA} !{FASTA_HEADER_DATA}
    python !{SCRIPT_DIR}/build_WGS_species_tree.py !{input} !{GENOME_DIR}
    echo "building gene trees"
    python !{SCRIPT_DIR}/generate_gene_trees.py !{GENE_FAMILIES} !{ALL_NUC_FASTA}
    echo "running markophylo"
    Rscript !{SCRIPT_DIR}/markophylo_indel_estiamtes.R !{SPECIES_TREE_FILE} !{PA_MATRIX} !{ALL_DATA} 2>| markophylo.errs
    echo "parsing markophylo"
    python !{SCRIPT_DIR}/parse_markophylo_output.py ./
    echo "building  network"
    Rscript !{SCRIPT_DIR}/network_builder.R ./
    """
}

