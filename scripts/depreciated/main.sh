#!/bin/bash
start="$(date +%s)"
echo "start time: $(date)"
scriptdir="/home/sid/thesis_SidReed/scripts"
crispr_annotation_info="/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json"
genus="$1" #genus name, capitalize first lettter
echo "fetching/formatting sequences for $genus"
mkdir -p "$genus"
cd "$genus"
python "$scriptdir"/get_genus_fastas.py "$genus" "$crispr_annotation_info"
#nucleotide
python "$scriptdir"/shortenFastaheaders.py nucleotide_fastas
cat nucleotide_fastas/*.fna >| all_nucleotides.fna
#protein
python "$scriptdir"/shortenFastaheaders.py protein_fastas
cat protein_fastas/*.faa >| all_proteins.faa
echo "single link clustering to gene families"
mkdir -p ./mmseq_clustering
cd ./mmseq_clustering
"$scriptdir"/mmseqWrapper.sh ../all_proteins.faa "$genus"
cd ..
python "$scriptdir"/parse_mmseq_clusters.py ./mmseq_clustering/clusters_"$genus".tsv ./all_proteins.faa
echo "creating PA matrix"
python "$scriptdir"/createPAMatrix.py gene_families.json protein_fastas/
echo "building species trees"
python "$scriptdir"/build_16S_species_tree_from_families.py all_nucleotides.fna fasta_headers_info.json
echo "running markophylo"
Rscript "$scriptdir"/markophylo_indel_estiamtes.R "species_tree_files/species_tree/species_tree.con.tre" binary_pa_matrix.csv "$crispr_annotation_info" 2>| markophylo.errs
echo "parsing markophylo"
python "$scriptdir"/parse_markophylo_output.py ./
echo "building gene trees"
python "$scriptdir"/generate_gene_trees.py gene_families.json all_nucleotides.fna
echo "building  network"
Rscript "$scriptdir"/network_builder.R ./
echo "generating network report"
python "$scriptdir"/network_analysis.py ./
end="$(date +%s)"
echo "end time: $(date)"
hrs="$(echo "($end-$start)/3600" | bc)"
mins="$(echo "($end-$start-$hrs*3600)/60" | bc)"
echo "runtime for $genus: $hrs hours $mins minutes"

