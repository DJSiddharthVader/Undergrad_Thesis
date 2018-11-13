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
    genelist = []
    for seq in SeqIO.parse(allgenesfasta,'fasta'):
        genelist.append(seq.id)
    return list(set(genelist))

def filterblasttable(data,coverage,pident,evalue):
#this should already be filtered by diamond, but just in case
    #allgenes=tuple(set(data['qseqid'].values.tolist() + data['sseqid'].values.tolist()))
    newdata = data[data['pident'] > pident]
    newdata = newdata[(newdata['qlen']/newdata['slen'] > coverage)]
    newdata = newdata[newdata['qseqid'] != newdata['sseqid']]
    newdata = newdata[newdata['evalue'] < evalue]
    return newdata#,allgenes

def fasterbuilddmat(data,genes):
    geneidxs = {gene:i for i,gene in enumerate(genes)}
    dmat = np.ones((len(genes),len(genes)))
    for n,row in tqdm(data.iterrows(),total=len(data)):
        i = geneidxs[row['qseqid']]
        j = geneidxs[row['sseqid']]
        dmat[i][j] = 0
        dmat[j][i] = 0
    for i in range(len(genes)):
        dmat[i][i] = 0 #self hits suppressed, defualt dist is 1, set diagonal to 0 manual
    idxgenes = {v:k for k,v in geneidxs.items()}
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
    singlerecords = []
    with open(fastaname) as fasta:
        for record in SeqIO.parse(fasta,'fasta'):
            if record.id in genelist:
                singlerecords.append(record)
    foutname = 'singletons_{}.faa'.format(outname)
    with open(foutname,'w') as singleout:
        SeqIO.write(singlerecords,foutname,'fasta')
    return foutname

def main2(ogdata,allgenesfasta,outname):
    genes = getgenelist(allgenesfasta)
    data= filterblasttable(ogdata,0.85,85,0.05)
    print('building distance matrix...')
    dmat,idxgenes = fasterbuilddmat(data,genes)
    print('{} unique genes'.format(len(list(idxgenes))))
    print('clustering...')
    fclust = linkage(dmat)
    families = format_to_dict2(fclust,idxgenes)
    families,singletons = filtersingletons(families)
    print('fetching singleton families...')
    foutname = singlegenelisttofasta(singletons,allgenesfasta,outname)
    print('{} singletons genes'.format(len(singletons)))
    return families,foutname

if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1],delimiter='\t') #blast/diamond table
    families,foutname = main2(data,sys.argv[2],sys.argv[3])
    json.dump(families,open(sys.argv[3] + '.json','w'))
    #json.dump(singletons,open('singletons_' + sys.argv[3],'w'))
    print('data output to file {}.json'.format(sys.argv[3]))
    print('singlton families in {}'.format(foutname))

#old stuff
#def format_to_dict(fclust,genes):
#    families = defaultdict(list)
#    for n,i in enumerate(fclust):
#        families['fam' + str(i)].append(genes[n])
#    return families
#def builddmat(data):
#    genes = tuple(set(data['qseqid'].values.tolist() + data['sseqid'].values.tolist()))
#    dmat = np.ones((len(genes),len(genes)))
#    for i,gi in enumerate(tqdm(genes)):
#        for j,gj in enumerate(tqdm(genes)):
#            dmat[i][j] = bindist(gi,gj,data) #every match has dist 0
#    return dmat,genes
#singlefams = defaultdict(list)
#for i,r in tqdm(data.iterrows(),total=len(data)):
#    singlefams[r['qseqid']].append(r['sseqid'])
#print(len(list(singlefams.keys())) == len(list(set(data['qseqid']))))
#print(set([len(x) for x in singlefams.values()]))
#singlefams = [[x] + f for x,f in zip(singlefams.keys(),singlefams.values())] #add index member to each family now list of lists, each list is all genes that have a blast hit (blast already filtered when run
#merged = True
#while merged:
#    merged = False
#    for idx,fam in enumerate(tqdm(singlefams)):
#        for i2,fam2 in enumerate(singlefams):
#            if set(fam).intersection(fam2) != None:
#                singlefams[idx] = list(set(fam).union(fam2))
#                merged = True
#json.dump(singlefams,open(sys.argv[2],'w'))
#def main(data):
#    dmat,genes = builddmat(data)
#    fclust = linkage(dmat)
#    families = format_to_dict(fclust,genes)
#    return families
