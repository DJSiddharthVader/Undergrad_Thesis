import sys
import json
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
from tqdm import tqdm
import fastcluster as fc
from functools import partial
from collections import defaultdict
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from multiprocessing.dummy import Pool as ThreadPool

#arg 1 is blast table
#arg 2 is concated all genes protein fasta
#arg 3 is genus name

def filterblasttable(data,coverage,pident,evalue):
#this should already be filtered by diamond, but just in case
    #allgenes=tuple(set(data['qseqid'].values.tolist() + data['sseqid'].values.tolist()))
    newdata = data[data['pident'] > pident]
    newdata = newdata[(newdata['qlen']/newdata['slen'] > coverage)]
    newdata = newdata[newdata['qseqid'] != newdata['sseqid']]
    newdata = newdata[newdata['evalue'] < evalue]
    return newdata

def fasterbuilddmat(data,genes):
    geneidxs = {gene:i for i,gene in enumerate(genes)}
    idxgenes = {v:k for k,v in geneidxs.items()}#same as geneidxs, reversed
    dmat = np.ones((len(genes),len(genes))) #set everything to have distance of 1 (max difference)
    for n,row in tqdm(data.iterrows(),total=len(data)):
        i = geneidxs[row['qseqid']]
        j = geneidxs[row['sseqid']]
        dmat[i][j] = 0
        dmat[j][i] = 0
    for i in range(len(genes)):
        dmat[i][i] = 0 #self hits suppressed, defualt dist is 1, set diagonal to 0 manual
    return dmat,idxgenes

def fastlinkage(dmat):
    return sch.fcluster(fc.linkage(squareform(dmat),method='single'),0.01)

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
    fclust = fastlinkage(dmat)
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

#depreciated
def fsize(fname):
    p = subprocess.Popen(['wc','-l',fname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def buildMatChunk(chunk,genes,coverage,pident,evalue):
    chunk = filterblasttable(chunk,converage,pident,evalue)
    newfile = open('tmp_dmat.npy','wb')
    for n,row in tqdm(chunk.iterrows(),total=chunk.shape[0],leave=False):
            do thing
    return dmat

def buildDistMatParallel(genes,table,coverage,pident,evalue,chunksize=100000):
    chunks = pd.read_csv(table,sep='\t',chunksize=chunksize)
    dmatfnc = partial(buildMatChunk,genes=genes,
                                    coverage=coverage,
                                    pident=pident,
                                    evalue=evalue)
    total = fsize(table)/chunksize
    outfiles = list(tqdm(pool.imap(dmatfnc,chunks),total=total))
    return outfiles

def getgenelist(allgenesfasta):
    return list(set([seq.id for seq in SeqIO.parse(allgenesfasta,'fasta')]))

#benchmarking
def timer(fnc):
    import time
    def wrapper(*args,**kwargs):
        s = time.time()
        result = fnc(*args,**kwargs)
        e = time.time()
        print('took {} time'.format(e -s))
        return result
    return wrapper

@timer
def prep(blastable,allfasta):
    genes = getgenelist(allfasta)
    data = filterblasttable(blastable,0.85,85,0.05)
    dmat,idxgenes = fasterbuilddmat(data,genes)
    square = squareform(dmat) #needed for linkage methdos
    return square

@timer
def scipyLinkage(dmat):
    return sch.fcluster(sch.single(dmat),0.01)

@timer
def fcLinkage(dmat):
    return sch.fcluster(fc.linkage(dmat,method='single'),0.01)

def compare(blastable,allfasta):
    square = prep(blastable,allfasta)
    spy = scipyLinkage(square)
    fcl = fcLinkage(square)
    print(np.array_equal(spy,fcl))
    return None

def linkage(dmat):
    square = squareform(dmat) #needed for linkage methdos
    linkmat = sch.single(square)
    return sch.fcluster(linkmat,0.0001)
