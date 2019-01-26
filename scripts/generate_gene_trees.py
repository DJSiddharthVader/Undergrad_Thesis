import os
import sys
import json
import math
import glob
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from functools import partial
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO,AlignIO,Alphabet
from multiprocessing.dummy import Pool as ThreadPool

def getNumTaxa():
    return len(glob.glob('./protein_fastas/*.faa'))

def pickGeneFamiliesForTrees(genefamilies,minsize):
    #list of families that have exactly maxtax members and contain genes from every taxa being considered
    maxtax = getNumTaxa()
    if type(minsize) == float and (minsize > 0 and minsize < 1):
        minsize = int(math.floor(minsize*maxtax))
    print('minimum organisms gene family must be present in: {}'.format(minsize))
    return [fam for fam,members in tqdm(genefamilies.items(),total=len(genefamilies),desc='picking') if (len(set([x.split(':')[0] for x in members])) >= minsize)]

def extractSeqsForTree(headersNameTuple,nucFasta,outdir):
    genefam,headerlist = headersNameTuple
    seqRecordlist = []
    seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist] #put sequences in list
    fastafilename = '{}/fastas/{}.fna'.format(outdir,genefam)
    return fastafilename, seqRecordlist

def fixNames(seqlist):
    fixedseqs = []
    for sequence in seqlist:
        fixedseqs.append(SeqRecord(sequence.seq,
                                   name=sequence.name.split(':')[0],
                                   id=sequence.id.split(':')[0],
#                                   name=sequence.name.replace(':','__'),
#                                   id=sequence.id.replace(':','__'),
                                   description=''))
    return fixedseqs

def extractSeqsParallel(headerTuples,outdir,nucFasta):
    pool = ThreadPool(processes)
    extractfnc = partial(extractSeqsForTree,outdir=outdir,nucFasta=nucFasta)
    fastas = list(tqdm(pool.imap(extractfnc,headerTuples),total=len(headerTuples),desc='getseqs'))
    return fastas

def writeFastas(nameseqiter):
    for name,seqs in nameseqiter.items():
        with open(name,'w') as outfile:
            SeqIO.write(seqs,outfile,'fasta')#write all seqs to a fasta file
    return None

def alignGeneFamilies(fastafile,outdir):
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    fastabase = fastafile.split('.')[0]
    fastaname = '{}/fastas/{}'.format(outdir,fastafile)
    alnfile = '{}/tmp_{}.aln'.format(outdir,fastabase)
    align = '{} {} > {}'.format(mafftbase,fastaname,alnfile)
    os.system(align)
    AlignIO.convert(alnfile,'fasta','{}/nexus/{}.nex'.format(outdir,fastabase),'nexus',alphabet=Alphabet.generic_dna)
    os.remove(alnfile)
    return None

def parallelAlignFastas(outdir,processes):
    pool = ThreadPool(processes)
    fastas = os.listdir('{}/fastas'.format(outdir))
    alnfnc = partial(alignGeneFamilies,outdir=outdir)
    alnfiles = list(tqdm(pool.imap(alnfnc,fastas),total=len(fastas),desc='aligning'))
    return None

def buildTree(nexfile,outdir):
    nexbase = nexfile.split('.')[0]
    mbf = """set autoclose=yes nowarn=yes
execute {}/nexus/{}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file={}
sump burnin=250
sumt burnin=250
quit""".format(outdir,nexfile,nexbase)
    mbscriptname = '{}_mrbayes_script.txt'.format(nexbase)
    open(mbscriptname,'w').write(mbf)
    logfilename = '{}_mrbayes_species_tree_log.txt'.format(nexbase)
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    treedir = '{}/trees/{}'.format(outdir,nexbase)
    os.system('mkdir -p {}'.format(treedir))
    for fp in glob.glob('./{}*'.format(nexbase)):
        os.rename(fp,os.path.join(os.getcwd(),treedir,fp))
    return None

def parallelBuildTrees(outdir,processes):
    pool = ThreadPool(processes)
    alnfiles = os.listdir('{}/nexus'.format(outdir))
    treefnc = partial(buildTree,outdir=outdir)
    treefiles = list(tqdm(pool.imap(treefnc,alnfiles),total=len(alnfiles),desc='trees'))
    #test = tqdm(pool.map(treefnc,alnfiles))
    return None

def main(familys,nucFasta,genusname,minsize,processes):
    #set up out directories
    outdir = 'gene_tree_files'
    os.system('mkdir -p {}'.format(outdir))
    os.system('mkdir -p {}/fastas'.format(outdir))
    os.system('mkdir -p {}/nexus'.format(outdir))
    os.system('mkdir -p {}/trees'.format(outdir))
    #pick families
    maxtax = getNumTaxa()
    family_idxs = pickGeneFamiliesForTrees(familys,minsize)
    print('Using {} of {} total families'.format(len(family_idxs),len(familys)))
    headerlists = {famidx:familys[famidx] for famidx in family_idxs}
    headerTuples = [(famidx,familys[famidx]) for famidx in family_idxs]
    #write families to fasta files
    seqs_to_write = extractSeqsParallel(headerTuples,outdir,nucFasta)
    seqs_to_write = {name:fixNames(seqs) for (name,seqs) in seqs_to_write}
    writeFastas(seqs_to_write)
    #align gene family fastas and convert to nexus
    parallelAlignFastas(outdir,processes)
    #build gene trees with mrbayes from alingments
    parallelBuildTrees(outdir,processes)
    return None


if __name__ == '__main__':
    famidx = json.load(open(sys.argv[1]))
    nucFasta = sys.argv[2]
    genusname = sys.argv[3]
    minsize = 0.4
    processes = 48
    main(famidx,nucFasta,genusname,minsize,processes)

#DEPRECIATED
def pickGeneFamiliesForTrees(pamat,columnindex,minsize):
    gene_family_list = []
    if type(minsize) == float and (minsize > 0 and minsize < 1):
        minsize = int(math.floor(minsize*pamat.shape[0]))
    print('minimum organisms gene family must be present in: {}'.format(minsize))
    binmat = matrix_to_binary(pamat)
    sizes = np.apply_along_axis(sum,0,binmat)
    family_idxs = [str(i) for i,size in enumerate(sizes) if size >= minsize]
    return family_idxs
