import os
import re
import json
import pandas as pd
from functools import reduce

def loaddf(dfpath):
    return pd.DataFrame(json.load(open(dfpath)))

def prepCrispr(path):
    df = loaddf(path)
    df['align'] = df['GCA'].apply(lambda x: '_'.join(x.split('_')[1:]))
    df['Genus'] = df['Name'].apply(lambda x: re.sub('[^a-zA-Z]','',x.split(' ')[0]))
    df['isCRISPR'] =df['System Types (CRISPRone)'].apply(lambda x: x not in ['NA','No System'])
    return(df)

def prepNCBI(path):
    df = loaddf(path)
    df['align'] = df['GCA'].apply(lambda x: '_'.join(x.split('_')[1:]))
    df['acctrim'] = df['Accession Number'].apply(lambda x: x.split('.')[0])
    return df

def prepDownloads():
    protpath = '/home/sid/thesis_SidReed/data/prot_cds_genomic'
    nucsuffix = '_cds_from_genomic.fna.gz'
    nucpath = '/home/sid/thesis_SidReed/data/nuc_cds_genomic'
    protsuffix = '_translated_cds.faa.gz'
    gcanames = ['_'.join(x.split('_')[:3]) for x in os.listdir(protpath)]
    df = pd.DataFrame({'gcanames':gcanames})
    df['protPath'] = df['gcanames'].apply(lambda x:'{}/{}{}'.format(protpath,x,protsuffix))
    df['nucPath'] = df['gcanames'].apply(lambda x:'{}/{}{}'.format(nucpath,x,nucsuffix))
    df['align'] = df['gcanames'].apply(lambda x: '_'.join(x.split('_')[1:3]))
    return(df)

def main(crisprpath,ncbipath,outpath):
    cdf = prepCrispr(crisprpath)
    ncbi = prepNCBI(ncbipath)
    ddf = prepDownloads()
    dfs = [cdf,ncbi,ddf]
    cn = pd.merge(ncbi,ddf,on='align',how='left')
    print(cn.shape)
    fdf = pd.merge(cn,cdf,on='align',how='inner')
    print(fdf.shape)
    #fdf = reduce(lambda left,right: pd.merge(left,right,on='align',how='outer'), dfs)
    fdf.to_json(outpath)
    return None

if __name__ == '__main__':
    crisprpath = '/home/sid/thesis_SidReed/data/allCRISPRAnnotationData/CRISPRone_files/mp_CRISPRoneAnnotations_01_11_18SR.json'
    ncbipath = '/home/sid/thesis_SidReed/data/updated_NCBI_info_16_03_19.json'
    outpath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'
    main(crisprpath,ncbipath,outpath)
