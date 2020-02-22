
//Globals
DATA_DIR="/home/sid/thesis_SidReed/data/"
SCRIPT_DIR="/home/sid/thesis_SidReed/scripts"
PROCESSES=8

//Input Files
GENOMIC_DIR="${DATA_DIR}/genomic_data"
ALL_DATA="${DATA_DIR}/all_data.tsv"
NUC_DIR="nucleotide"
PROT_DIR="protein"
GENOME_DIR="genome"

//Output Files/Dirs
MMSEQ_DIR="mmseq_clustering"
FASTA_HEADER_DATA="fasta_headers_info.json"
GENE_FAMILIES="gene_families.json"
ALL_NUC_FASTA="all_nucleotides.fna"
ALL_PROT_FASTA="all_proteins.faa"
PA_MATRIX="binary_pa_matrix.csv"
SPECIES_TREE_FILE="species_tree_WGS/WGS_species_tree.con.tree"


setUpFiles = {
    // $input is a genus name
    // happens in $DATA_DIR
    exec """
        mkdir -p "$input" && cd "$input"
        python "$SCRIPT_DIR"/get_genus_fastas.py "$ALL_DATA" ./ "$input" "$GENOMIC_DIR"
        python "$SCRIPT_DIR"/shortenFastaheaders.py "$NUC_DIR"
        cat "$NUC_DIR"/*.fna >| "$ALL_NUC_FASTA"
        python "$SCRIPT_DIR"/shortenFastaheaders.py "$PROT_DIR"
        cat "$PROT_DIR"/*.faa >| "$ALL_PROT_FASTA"
    """
}
clusterGenes = {
    exec """
        echo "single link clustering to gene families"
        "$SCRIPT_DIR"/mmseqWrapper.sh "$ALL_PROT_FASTA" "$input" "$MMSEQ_DIR"
        python "$SCRIPT_DIR"/parse_mmseq_clusters.py ./mmseq_clustering/clusters_"$input".tsv "$ALL_PROT_FASTA"
        echo "creating P/A matrix"
        python "$SCRIPT_DIR"/createPAMatrix.py "$GENE_FAMILIES" "$PROT_DIR"
    """
}
buildTrees = {
    exec """
        echo "building species trees"
        // python "$SCRIPT_DIR"/build_16S_species_tree_from_families.py "$ALL_NUC_FASTA" "$FASTA_HEADER_DATA"
        python "$SCRIPT_DIR"/build_WGS_species_tree.py "$input" "$GENOME_DIR"
        echo "building gene trees"
        python "$SCRIPT_DIR"/generate_gene_trees.py "$GENE_FAMILIES" "$ALL_NUC_FASTA"
    """
}
generateResults = {
    exec """
        echo "running markophylo"
        Rscript "$SCRIPT_DIR"/markophylo_indel_estiamtes.R "$SPECIES_TREE_FILE" "$PA_MATRIX" "$ALL_DATA" 2>| markophylo.errs
        echo "parsing markophylo"
        python "$SCRIPT_DIR"/parse_markophylo_output.py ./
        echo "building  network"
        Rscript "$SCRIPT_DIR"/network_builder.R ./
    """
}

run = {
    setUpFiles + clusterGenes + buildTrees //+ generateResults
}
