import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from collections import defaultdict
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import single,fcluster

#arg 1 is blast table
#arg 2 is concated all genes fasta
#arg 3 is output name (no extension)

def getgenelist(allgenesfasta):
    return list(set([seq.id for seq in SeqIO.parse(allgenesfasta,'fasta')]))

def filterblasttable(data,coverage,pident,evalue):
#this should already be filtered by diamond, but just in case
    #allgenes=tuple(set(data['qseqid'].values.tolist() + data['sseqid'].values.tolist()))
    newdata = data[data['pident'] > pident]
    newdata = newdata[(newdata['qlen']/newdata['slen'] > coverage)]
    newdata = newdata[newdata['qseqid'] != newdata['sseqid']]
    newdata = newdata[newdata['evalue'] < evalue]
    return newdata

def fasterbuilddmat(data,genes):
    geneidxs = {gene:i for i,gene in enumerate(genes)} #create row indices for each gene
    idxgenes = {v:k for k,v in geneidxs.items()}#same as geneidxs but reverse
    dmat = np.ones((len(genes),len(genes))) #set everything to have distance of 1 (max difference)
    for n,row in tqdm(data.iterrows(),total=len(data)):
        i = geneidxs[row['qseqid']]
        j = geneidxs[row['sseqid']]
        dmat[i][j] = 0
        dmat[j][i] = 0
    for i in range(len(genes)):
        dmat[i][i] = 0 #self hits suppressed, defualt dist is 1, set diagonal to 0 manual
    return dmat,idxgenes

def linkage(dmat):
    square = squareform(dmat) #needed for linkage methdos
    linkmat = single(square)
    return fcluster(linkmat,0.0001)

def format_to_dict2(fclust,idxgenes):
    families = defaultdict(list)
    for n,i in enumerate(fclust):
        families['fam' + str(i)].append(idxgenes[n])
    return families

def filtersingletons(families):
    singletons = []
    newfams = {}
    for i,fam in families.items():
        if len(fam) == 1:
            singletons.append(fam[0])
        else:
            newfams[i] = fam
    return newfams,singletons

def singlegenelisttofasta(genelist,fastaname,outname):
    with SeqIO.parse(open(fastaname),'fasta') as fasta:
        singlerecords = [seq for seq in fasta if seq.id in genelist]
    foutname = 'singletons_{}.faa'.format(outname)
    with open(foutname,'w') as singleout:
        SeqIO.write(singlerecords,foutname,'fasta')
    return foutname

def main(ogdata,allgenesfasta,outname):
    genes = getgenelist(allgenesfasta)
    data = filterblasttable(ogdata,0.85,85,0.05)
    print('building distance matrix...')
    dmat,idxgenes = fasterbuilddmat(data,genes)
    print('{} unique genes'.format(len(list(idxgenes))))
    print('clustering...')
    fclust = linkage(dmat)
    families = format_to_dict2(fclust,idxgenes)
    print('fetching singleton families...')
    families,singletons = filtersingletons(families)
    foutname = singlegenelisttofasta(singletons,allgenesfasta,outname)
    print('{} singletons genes'.format(len(singletons)))
    return families,foutname

if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1],delimiter='\t') #blast/diamond table
    families,foutname = main(data,sys.argv[2],sys.argv[3])
    json.dump(families,open('gene_families_{}.json'.format(sys.argv[3]),'w'))
    #json.dump(singletons,open('singletons_' + sys.argv[3],'w'))
    print('data output to file {}.json'.format(sys.argv[3]))
    print('singlton families in {}'.format(foutname))

