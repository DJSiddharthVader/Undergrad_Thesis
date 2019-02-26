#!/bin/bash
scriptdir="/home/sid/thesis_SidReed/scripts"
fps="$1" #file with abs paths to all genomes being used
genusname="$2" #genus name
echo "preping sequences"
"$scriptdir"/prepSeqs.sh "$fps" "$genusname"
echo "single link clustering to gene families"
cd "$genusname"
mkdir mmseq_clustering
cd mmseq_clustering
"$scriptdir"/mmseqWrapper.sh ../all_proteins.faa "$genusname"
cd -
python "$scriptdir"/parse_mmseq_clusters.py ./mmseq_clustering/clusters_"$genusname".tsv ./all_proteins.faa
echo "creating PA matrix"
python "$scriptdir"/createPAMatrix.py gene_families.json protein_fastas/ "$genusname"
echo "building species trees"
python "$scriptdir"/build_species_tree.py gene_families.json all_nucleotides.fna "$genusname"
echo "building gene trees"
python "$scriptdir"/generate_gene_trees.py gene_families.json all_nucleotides.fna "$genusname"
python "$scriptdir"/build_network_from_trees.py ../"$genusname"

#old main
#echo "preping sequences"
#"$scriptdir"/prepSeqs.sh "$fps" "$genusname"
#echo "all vs all alingment"
#"$scriptdir"/diamondAlign.sh "$genusname"
#cd "$genusname"
#echo "single link clustering to gene families"
#python "$scriptdir"/single_link_clustering_blast_hits.py allvsallproteins.dmnd.out all_proteins.faa "$genusname"
#echo "creating PA matrix"
#python "$scriptdir"/createPAMatrix.py gene_families.json protein_fastas/ "$genusname"
#echo "building species trees"
#python "$scriptdir"/build_species_tree.py pa_matrix_"$genusname".npy column_indexes_families_"$genusname".json gene_families_"$genusname".json all_nucleotides.fna "$genusname"
#echo "building gene trees"
#python "$scriptdir"/generate_gene_trees.py pa_matrix_"$genusname".npy column_indexes_families_"$genusname".json gene_families_"$genusname".json all_nucleotides.fna "$genusname"
#python "$scriptdir"/build_network_from_trees.py ../"$genusname"
