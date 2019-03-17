#!/bin/python3
import os
import sys
import json
import pandas as pd

#df_path = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'

def mkdirs(destdir):
    nucdir = "{}/nucleotide_fastas".format(destdir)
    os.system("mkdir -p {}".format(nucdir))
    protdir = "{}/protein_fastas".format(destdir)
    os.system("mkdir -p {}".format(protdir))
    return(nucdir,protdir)

def copyfiles(df_path,genus,nucdir,protdir):
    nuc_suffix = "_cds_from_genomic.fna.gz"
    prot_suffix = "_translated_cds.faa.gz"
    meta_data_df = pd.DataFrame(json.load(open(df_path)))
    genus_df = meta_data_df[meta_data_df['Genus'] == genus]
    for i,row in genus_df.iterrows():
        cpnuc = "cp {} {}/{}{}".format(row['nucPath'],
                                       nucdir,
                                       row['Accession Number'],
                                       nuc_suffix)
        os.system(cpnuc)
        cpprot = "cp {} {}/{}{}".format(row['protPath'],
                                        protdir,
                                        row['Accession Number'],
                                        prot_suffix)
        os.system(cpprot)
    os.system("gzip -d {}/*.gz".format(nucdir))
    os.system("gzip -d {}/*.gz".format(protdir))

def main(genus,df_path,destdir):
    nucdir,protdir = mkdirs(destdir)
    copyfiles(df_path,genus,nucdir,protdir)
    return None

if __name__ == '__main__':
    genus = sys.argv[1]
    df_path = sys.argv[2]
    if len(sys.argv) < 3:
        destdir = sys.argv[3]
    else:
        destdir = os.getcwd()
    main(genus,df_path,destdir)



#DEPRECIATED
#suffixes and directories
#source_nuc_dir = "/home/sid/thesis_SidReed/data/nuc_cds_genomic"
#source_prot_dir = "/home/sid/thesis_SidReed/data/prot_cds_genomic"
