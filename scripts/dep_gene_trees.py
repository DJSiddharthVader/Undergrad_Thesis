import os
import sys
import json
import math
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.Nexus import Nexus
from Bio import SeqIO,AlignIO,Alphabet
from multiprocessing.dummy import Pool as ThreadPool

#pick genes families with col len([x >0 for x in a]) > len(col)*-.40, what coverage cutoff?
#get one or all members, from each family?
#build fasta
#align fasta
#build tree
#write tree to  file

def pickFamiliesForTree(pamat,columnindex,presence_cutoff=0.3):
    famlist = []
    coff = int(math.floor(presence_cutoff*pamat.shape[0])) #only get families present in cutoff% of organisms
    for  i,family in enumerate(pamat.T):
        if len([1 for x in family if x > 0]) > coff:
            famlist.append(columnindex[str(i)])
    return famlist

def extractSeqsForTree(headerlist,allfasta):
    seqRecordlist = []
    seqiter = SeqIO.parse(open(allfasta),'fasta') #open fasta file
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist] #put sequences in list
    return seqRecordlist

def alnseqs(fastapath,outname)
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    alnfile = 'tmp_{}.aln'.format(outname)
    align = '{} {} > {}'.format(mafftbase,fastapath,alnfile)
    os.system(align)

def buildTree(fastapath,outdir,outname,fam):
    mbf =\
"""set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file={}/{}_{}
quit""".format(alnfile,'gene_tree',outdir,outname,fam)
    mbscriptname = 'mbScript_{}.txt'.format(outname)
    open(mbscriptname,'w').write(mbf)
    logfilename = 'mb_gene_tree_log_{}.txt'.format(outname)
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    return logfilename

def main(pamat,colidx,famdict,allnucfasta,outname):
    genefams = pickFamiliesForTree(pamat,colidx)
    print('using {}/{} genes for trees'.format(len(genefams),pamat.shape[1]))
    fadir = 'gene_tree_fastas_{}'.format(outname)
    os.system('mkdir -p {}'.format(fadir))
    treedir = 'mb_gene_trees_{}'.format(outname)
    os.system('mkdir -p {}'.format(treedir))
    for fam in tqdm(genefams):
        headerlist = famdict[fam] #all genes in the family
        seqs = extractSeqsForTree(headerlist,allnucfasta) #get all sequences in the family
        fasta = '{}/{}_{}.fna'.format(fadir,outname,fam)
        with open(fasta,'w') as out:
            SeqIO.write(seqs,fasta,'fasta') #write all seqs to a fasta file
        treepath = buildTree(fasta,treedir,outname,fam)
    for fp in os.listdir(treedir):
        os.rename(fp,os.path.join(os.getcwd(),treedir,fp))
    return treedir

if __name__ == '__main__':
    pamat = np.load(sys.argv[1])
    columnindex = json.load(open(sys.argv[2]))
    familyindex = json.load(open(sys.argv[3]))
    nucFasta = sys.argv[4]
    outname = sys.argv[5]
    main(pamat,columnindex,familyindex,nucFasta,outname)

