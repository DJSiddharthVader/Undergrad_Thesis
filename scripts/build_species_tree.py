import os
import sys
import json
import glob
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.Nexus import Nexus
from functools import partial
from Bio import SeqIO,AlignIO,Alphabet
from multiprocessing.dummy import Pool as ThreadPool

#arg 1 is PA matrix (not binary)
#arg 2 is col families index dict
#arg 3 is actualy genefamilies dict
#arg 4 is concatenated fasta fle with all genes (nuc seqs)
#arg 5 is outname, no extension

#will extract all genes present in only 1 copy in every taxa and spit to a fasta file for alignment and tree creation

def matrix_to_binary(pamat):
    binarize = lambda x:np.where(x > 0, 1, 0)
    vecbin = np.vectorize(binarize)
    return vecbin(pamat)

def pickGeneFamilies(pamat,columnindex):
    binmat = matrix_to_binary(pamat)
    sizes = np.apply_along_axis(sum,0,binmat)
    family_idxs = [str(i) for i,size in enumerate(sizes) if size == binmat.shape[0]]
    print('using {} of {} genes for species tree'.format(len(family_idxs),binmat.shape[1]))
    return family_idxs

def extractSeqsForTree(headerlist,nucFasta,outdir,genefam):
    seqRecordlist = []
    seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist] #put sequences in list
    fastafilename = '{}/fastas/{}.fna'.format(outdir,genefam)
    return fastafilename, seqRecordlist

def fixNames(seqlist):
    fixedseqs = []
    for sequence in seqlist:
        fixedseqs.append(SeqRecord(sequence.seq, name=sequence.name.split(':')[0], id=sequence.id.split(':')[0]))
    return fixedseqs

def writeFastas(nameseqiter):
    for name,seqs in nameseqiter.items():
        with open(name,'w') as outfile:
            SeqIO.write(seqs,outfile,'fasta') #write all seqs to a fasta file
    return None

def concatNexAlns(nexDir,outname,same_taxa=True):#from https://biopython.org/wiki/Concatenate_nexus
    """Combine multiple nexus data matrices in one partitioned file.
    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this """
    filelist = [x for x in os.listdir(nexDir) if x.endswith('.nex')]
    nexi = [(os.path.join(nexDir,fname), Nexus.Nexus(os.path.join(nexDir,fname))) for fname in filelist]
    coutname = 'concat_stree_aln_{}.nex'.format(outname)
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open(coutname,'w'))
    return coutname

def buildMbSpecTree(nexaln,outname):
    mbf =\
"""set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file={}_{}
quit""".format(nexaln,'spectree',outname)
    mbscriptname = 'mbScript_{}.txt'.format(outname)
    open(mbscriptname,'w').write(mbf)
    logfilename = 'mb_spectree_log_{}.txt'.format(outname)
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    outdir = 'mb_species_tree_{}'.format(outname)
    os.system('mkdir -p {}'.format(outdir))
    for fp in glob.glob('./*spectree*'):
        os.rename(fp,os.path.join(os.getcwd(),outdir,fp))
    os.rename(mbscriptname,os.path.join(os.getcwd(),outdir,mbscriptname))
    return outdir

def main(pamat,columns,familys,nucFasta,genusname,genes):
    #set up out directories
    outdir = 'species_tree_files'
    os.system('mkdir -p {}'.format(outdir))
    os.system('mkdir -p {}/fastas'.format(outdir))
    os.system('mkdir -p {}/nexus'.format(outdir))
    os.system('mkdir -p {}/trees'.format(outdir))
    #pick families
    family_col_idxs = pickGeneFamiliesForTrees(pamat,columns,minsize)[:genes]
    family_idxs = [columns[col_idx] for col_idx in family_col_idxs]
    headerlists = {famidx:familys[famidx] for famidx in family_idxs}
    #write families to fasta files
    seqs_to_write = [extractSeqsForTree(heads,nucFasta,outdir,fam) for fam,heads in tqdm(headerlists.items(),total=len(list(headerlists.keys())),desc='getseqs')]
    seqs_to_write = {name:fixNames(seqs) for (name,seqs) in seqs_to_write}
    writeFastas(seqs_to_write)
    #align gene family fastas and convert to nexus
    alignGeneFamilies(outdir)
    #buiil gene trees with mrbayes from alingments
    buildTrees(outdir)
    return None

def main(pamat,columnindex,familyindex,nucFasta,outname,processes,fams):
    fastadir = 'species_tree_fastas_{}'.format(outname)
    fastadir,treegenes = writeAllSpeciesTreeFastas(pamat,columnindex,familyindex,nucFasta,fastadir,processes,fams)
    print('{} genes used for species tree'.format(treegenes))
    print('fasta files for species tree written to {}'.format(fastadir))
    nexDir = nexusAlnSTreeFastas(fastadir,fams)
    print('fasta alingments for species tree written to {}'.format(fastadir))
    print('nexus alingments for species tree written to {}'.format(nexDir))
    coutname = concatNexAlns(nexDir,outname,same_taxa=True)
    print('concat nexus aln file for mrbayes written to {}'.format(coutname))
    print('running MrBayes ...')
    toutname = buildMbSpecTree(coutname,outname)
    print('mrbayes species tree files in {}'.format(toutname))

if __name__ == '__main__':
    pamat = np.load(sys.argv[1])
    colidx = json.load(open(sys.argv[2]))
    famidx = json.load(open(sys.argv[3]))
    nucFasta = sys.argv[4]
    genusname = sys.argv[5]
    minsize = 0.3
    main(pamat,colidx,famidx,nucFasta,genusname,minsize)
