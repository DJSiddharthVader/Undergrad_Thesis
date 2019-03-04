#!/bin/python3
import os
import sys
import json
import pandas as pd

#suffixes and directories
nuc_suffix = "_cds_from_genomic.fna.gz"
source_nuc_dir = "/home/sid/thesis_SidReed/bacterialGBFFs/nuc_cds_genomic"
prot_suffix = "_translated_cds.faa.gz"
source_prot_dir = "/home/sid/thesis_SidReed/bacterialGBFFs/prot_cds_genomic"
df_path = '/home/sid/thesis_SidReed/pop_annotation_data_frame.json'

#args
genus = sys.argv[1]
if len(sys.argv) < 2:
    destdir = sys.argv[2]
else:
    destdir = os.getcwd()

#creatdirs
nucdir = "{}/nucleotide_fastas".format(destdir)
os.system("mkdir -p {}".format(nucdir))
protdir = "{}/protein_fastas".format(destdir)
os.system("mkdir -p {}".format(protdir))

#Copy,Rename fasta files
meta_data_df = pd.DataFrame(json.load(open(df_path)))
genus_df = meta_data_df[meta_data_df['Genus'] == genus]
for i,row in genus_df.iterrows():
    filename = 'GCF_' + '_'.join(row['GCA'].split('_')[1:])
    #print(filename)
    cpnuc = "cp {}/{}{} {}/{}{}".format(source_nuc_dir,
                                        filename,
                                        nuc_suffix,
                                        nucdir,
                                        row['Accession Number'],
                                        nuc_suffix)
    os.system(cpnuc)
    cpprot = "cp {}/{}{} {}/{}{}".format(source_prot_dir,
                                        filename,
                                        prot_suffix,
                                        protdir,
                                        row['Accession Number'],
                                        prot_suffix)
    os.system(cpprot)
os.system("gzip -d {}/*.gz".format(nucdir))
os.system("gzip -d {}/*.gz".format(protdir))
