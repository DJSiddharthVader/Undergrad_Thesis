import os
import sys
import json
import tempfile
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.Nexus import Nexus
from functools import partial
from Bio import SeqIO,AlignIO,Alphabet
from multiprocessing.dummy import Pool as ThreadPool

#pick genes families with col len([x >0 for x in a]) > len(col)*-.40, what coverage cutoff?
#get one or all members, from each family?
#build fasta
#align fasta
#build tree
#write tree to  file

def pickFamiliesForTree(pamat,colidx,cutoff=0.4,howmay=50):
    famlist = []
    coff = cutoff*len(pamat.shape[0])
    for  i,family in enumerate(pamat.T):
        if len([1 for x in family if x > 0]):
            famlist.append(colidx[str(i)])
    return famlist

def extractSeqsForTree(famlist,famdict,allfasta,firstonly=True):
    seqRecordlist = []
    seqiter = SeqIO.parse(open(allfasta),'w')
    headerlist = [famdict[fam][0] for fam in famlist]
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist]
    return seqRecordlist

def buildTree(fastapath,treedir,outname):
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    align ='{} {} > tmp_{}.aln'.format(mafftbase,fastapath,outname)
    os.system(align)

def main(pamat,colidx,famdict,allnucfasta,outname)
    genetreefams = pickFamiliesForTree(pamat,colidx)
    fadir = '{}_gtreeFastas'.format(outname)
    os.system('mkdir -p {}'.format(fadir))
    trdir = '{}_gtreeFastas'.format(outname)
    os.system('mkdir -p {}'.format(trdir))
    for t in list(range(trees):
        families = np.random.choice(genetreefams,50,replace=True)
        seqs = extractSeqsForTree(families,famdict,allnucfasta)
        fasta = '{}/{}_set{}.fna'.format(fadir,outname,t)
        with open(fasta,'w') as out:
            SeqIO.write(seqs,fasta,'fasta')
        treepath = buildTree(fasta,trdir,outname)
    return treedir

if __name__ == '__main__':
    pamat = np.load(sys.argv[1])
    columnindex = json.load(open(sys.argv[2]))
    familyindex = json.load(open(sys.argv[3]))
    nucFasta = sys.argv[4]
    outname = sys.argv[5]
    processes = 8
    fams = 50
    main(pamat,columnindex,familyindex,nucFasta,outname,processes,fams)
