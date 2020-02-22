import os
import re
import sys
import json
import pandas as pd

#Constants
columns = ['nuc_ftpurl','prot_ftpurl','genome_ftpurl']
suffixes = ['_cds_from_genomic.fna.gz', '_translated_cds.faa.gz','_genomic.fna.gz']
filenames = ['nucleotide','protein','genome',]


def loaddf(dfpath):
    return pd.DataFrame(json.load(open(dfpath)))

def prepCrispr(path):
    df = loaddf(path)
    df['align'] = df['GCA'].apply(lambda x: '_'.join(x.split('_')[1:])) #for mapping CRISPR One accessions to NCBI
    df['Genus'] = df['Name'].apply(lambda x: re.sub('[^a-zA-Z]','',x.split(' ')[0]))
    df['isCRISPR'] =df['System Types (CRISPRone)'].apply(lambda x: x not in ['NA','No System'])
    return(df)

def prepNCBI(ncbi_data,prok_data):
    df = loaddf(ncbi_data)
    pdf = pd.read_csv(prok_data)
    df = pd.merge(pdf,df,on='BioSample',how='inner')
    df['align'] = df['GCA'].apply(lambda x: '_'.join(x.split('_')[1:])) #for mapping CRISPR One accessions to NCBI
    df['acctrim'] = df['Accession Number'].apply(lambda x: x.split('.')[0]) #remove version from accession
    return df

def combineAndFilter(crispr,ncbi):
    merged_df = pd.merge(ncbi,crispr,on='align',how='inner') #only keep entries with CRISPR and NCBI info
    merged_df = merged_df[merged_df['RefSeq FTP'].notnull()] #remove entries without ftp links
    merged_df = merged_df[merged_df['Assembly level'] == 'Complete Genome'] #remove entries without ftp links
    return merged_df

def formatURL(url,suffix):
    filename = url.split('/')[-1]
    fileurl = '{}/{}{}'.format(url,filename,suffix)
    return fileurl

def writeURLs(df,filelist_dir):
    for column,suffix,fname in zip(columns,suffixes,filenames):
        df[column] = df['RefSeq FTP'].apply(lambda link: formatURL(link,suffix))
        outfile = os.path.join(filelist_dir,'{}_FTPURLs.txt'.format(fname))
        with open(outfile,'w') as filelist:
            filelist.write('\n'.join(df[column]))
    return df

def main(crispr_data,ncbi_data,prok_data,filelist_dir,outfile,min_strains_per_genera):
    crispr = prepCrispr(crispr_data)
    ncbi = prepNCBI(ncbi_data,prok_data)
    final = combineAndFilter(crispr,ncbi)
    final = writeURLs(final,filelist_dir)
    final['local_fname'] = final[columns[1]].apply(lambda x: os.path.basename(x).strip(suffixes[1]))
    genera_file = os.path.join(filelist_dir,'genera_greater_{}.txt'.format(min_strains_per_genera))
    with open(genera_file,'w') as filelist:
        genera = final['Genus'].value_counts()
        genera = list(genera[genera >= min_strains_per_genera].index)
        filelist.write('\n'.join(genera))
    #downloaded = prepDownloads(nuc_dir,prot_dir)
    #ncbi_crispr_df = combineAndFilter(crispr,ncbi,downloaded)
    #final = pd.merge(ncbi_crispr_df,downloaded,on='align',how='left')
    final.to_csv(outfile,sep='\t')
    print(crispr.shape)
    print(ncbi.shape)
    print(final.shape)
    return None


if __name__ == '__main__':
    crispr_data = sys.argv[1]
    ncbi_data = sys.argv[2]
    prok_data = sys.argv[3]
    filelist_dir= sys.argv[4]
    outfile = sys.argv[5]
    min_strains_per_genera = 5
    main(crispr_data,ncbi_data,prok_data,filelist_dir,outfile,min_strains_per_genera)

#DEPRECIATED
#def prepDownloads(nuc_dir,prot_dir):
#    gcanames = ['_'.join(x.split('_')[:3]) for x in os.listdir(prot_dir)]
#    df = pd.DataFrame({'gcanames':gcanames})
#    df['nuc_localfile'] = df['gcanames'].apply(lambda gca:'{}/{}{}'.format(nuc_dir,gca,nuc_suffix)) #add local filepaths for nucleotide files
#    df['prot_localfile'] = df['gcanames'].apply(lambda gca:'{}/{}{}'.format(prot_dir,gca,prot_suffix)) #add local filepaths for protein files
#    df['align'] = df['gcanames'].apply(lambda x: '_'.join(x.split('_')[1:3]))
#    return(df)
#
#    nucpath = '/home/sid/thesis_SidReed/data/nuc_cds_genomic'
#    prot_dir = '/home/sid/thesis_SidReed/data/prot_cds_genomic'
#    crisprpath = '/home/sid/thesis_SidReed/data/allCRISPRAnnotationData/CRISPRone_files/mp_CRISPRoneAnnotations_01_11_18SR.json'
#    ncbipath = '/home/sid/thesis_SidReed/data/updated_NCBI_info_16_03_19.json'
#    outpath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'
    #main(crisprpath,ncbipath,outpath)

