#!/bin/bash
start="$(date +%s)"
echo "start time: $start"
scriptdir="/home/sid/thesis_SidReed/scripts"
crispr_annotation_info="/home/sid/thesis_SidReed/data/pop_annotation_data_frame.json"
genus="$1" #genus name, capitalize first lettter
echo "fetching/formatting sequences"
mkdir -p "$genus"
cd "$genus"
python "$scriptdir"/get_genus_fastas.py "$genus"
#nucleotide
python "$scriptdir"/shortenFastaheaders.py nucleotide_fastas "$genus"
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
echo "building gene trees"
python "$scriptdir"/generate_gene_trees.py gene_families.json all_nucleotides.fna
echo "building  network"
Rscript "$scriptdir"/network_builder.R ./
end="$(date +%s)"
echo "end time: $end"
runtime=$((end-start))
echo "runtime: $runtime"

#DEPRECIATED
#old main
#echo "preping sequences"
#"$scriptdir"/prepSeqs.sh "$fps" "$genus"
#echo "all vs all alingment"
#"$scriptdir"/diamondAlign.sh "$genus"
#cd "$genus"
#echo "single link clustering to gene families"
#python "$scriptdir"/single_link_clustering_blast_hits.py allvsallproteins.dmnd.out all_proteins.faa "$genus"
#echo "creating PA matrix"
#python "$scriptdir"/createPAMatrix.py gene_families.json protein_fastas/ "$genus"
#echo "building species trees"
#python "$scriptdir"/build_species_tree.py pa_matrix_"$genus".npy column_indexes_families_"$genus".json gene_families_"$genus".json all_nucleotides.fna "$genus"
#echo "building gene trees"
#python "$scriptdir"/generate_gene_trees.py pa_matrix_"$genus".npy column_indexes_families_"$genus".json gene_families_"$genus".json all_nucleotides.fna "$genus"
#python "$scriptdir"/build_network_from_trees.py ../"$genus"
