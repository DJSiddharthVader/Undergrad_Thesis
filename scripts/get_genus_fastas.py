#!/bin/python3
import os
import sys
import pandas as pd

#initializes some directories and copies over the appropriate fasta files for a given genus
suffixes = ['_cds_from_genomic.fna.gz', '_translated_cds.faa.gz','_genomic.fna.gz']
datatypes = ['nucleotide','protein','genome']

def copyfiles(datatype,genus_df,genomic_dir,suffix):
    os.system("mkdir -p ./{}".format(datatype))#data/genus_data/$GENUS/nucleotide/
    for i,row in genus_df.iterrows():
        filename = '{}{}'.format(os.path.basename(row['local_fname']),#GCF_000005825.2_ASM582v2
                                 suffix) #_genomic.fna.gz
        cpcmd = "cp {}/{}/{} ./{}/{}".format(genomic_dir, #data/genomic_data
                                             datatype,    #nucleotide
                                             filename,    #GCF_000005825.2_ASM582v2_genomic.fna.gz
                                             datatype,    #nucleotide
                                             filename)    #GCF_000005825.2_ASM582v2_genomic.fna.gz
        #cp data/genomic_data/nucleotide/GCF_000005825.2_ASM582v2_genomic.fna.gz
        #   ./nucleotide/GCF_000005825.2_ASM582v2_genomic.fna.gz
        os.system(cpcmd)
    os.system("gzip -d ./{}/*.gz".format(datatype)) #unzip all files
    return None

def copyAllFiles(genus,df,genomic_dir):
    genus_df  = df[df['Genus'] == genus]
    for datatype,suffix in zip(datatypes,suffixes):
        copyfiles(datatype,genus_df,genomic_dir,suffix)
    return None

def main(df_path,genus,genomic_dir):
    df = pd.read_csv(df_path,sep='\t')
    copyAllFiles(genus,df,genomic_dir)
    return None


if __name__ == '__main__':
    df_path = sys.argv[1] #data/all_data.tsv
    genus = sys.argv[2] #the genus
    genomic_dir =  sys.argv[3] # data/genomic_data
    main(df_path,genus,genomic_dir)



#DEPRECIATED
#suffixes and directories
#source_nuc_dir = "/home/sid/thesis_SidReed/data/nuc_cds_genomic"
#source_prot_dir = "/home/sid/thesis_SidReed/data/prot_cds_genomic"
#initializes some directories and copies over the appropriate fasta files for a given genus
#nuc_suffix = "_cds_from_genomic.fna.gz"
#prot_suffix = "_translated_cds.faa.gz"
#
#def mkdirs(destdir):
#    #make directories for files to be downloaded
#    nucdir = "{}/nucleotide".format(destdir)
#    os.system("mkdir -p {}".format(nucdir))
#    protdir = "{}/protein".format(destdir)
#    os.system("mkdir -p {}".format(protdir))
#    return(nucdir,protdir)
#
#def copyfiles(df_path,genus,nucdir,protdir):
#    meta_data_df = pd.DataFrame(json.load(open(df_path)))
#    genus_df = meta_data_df[meta_data_df['Genus'] == genus]
#    for i,row in genus_df.iterrows():
#        cpnuc = "cp {} {}/{}{}".format(row['nucPath'], #copy,rename nuc fasta
#                                       nucdir,
#                                       row['Accession Number'],
#                                       nuc_suffix)
#        os.system(cpnuc)
#        cpprot = "cp {} {}/{}{}".format(row['protPath'], #copy,rename nuc fasta
#                                        protdir,
#                                        row['Accession Number'],
#                                        prot_suffix)
#        os.system(cpprot)
#    os.system("gzip -d {}/*.gz".format(nucdir)) #unzip all fasta files
#    os.system("gzip -d {}/*.gz".format(protdir))
#
#def main(genus,df_path,destdir):
#    nucdir,protdir = mkdirs(destdir)
#    copyfiles(df_path,genus,nucdir,protdir)
#    return None
#
#
#if __name__ == '__main__':
#    genus = sys.argv[1] #input
#    df_path = sys.argv[2] #input /home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json
#    if len(sys.argv) < 3:
#        destdir = sys.argv[3]
#    else:
#        destdir = os.getcwd() #output /home/sid/thesis_SidReed/data/genus_data/$GENUS
#    main(genus,df_path,destdir)
#
