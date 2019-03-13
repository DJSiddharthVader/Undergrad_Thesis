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

def extractSeqsForTree(headersNameTuple,nucFasta):
    genefam,headerlist = headersNameTuple
    seqRecordlist = []
    seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist] #put sequences in list
    fastafilename = 'gene_tree_files/fastas/{}.fna'.format(genefam)
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

def extractSeqsParallel(headerTuples,nucFasta):
    pool = ThreadPool(processes)
    extractfnc = partial(extractSeqsForTree,nucFasta=nucFasta)
    fastas = list(tqdm(pool.imap(extractfnc,headerTuples),total=len(headerTuples),desc='getseqs'))
    return fastas

def writeFastas(nameseqiter):
    for name,seqs in nameseqiter.items():
        with open(name,'w') as outfile:
            SeqIO.write(seqs,outfile,'fasta')#write all seqs to a fasta file
    return None

def alignGeneFamilies(fastafile):
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    fastabase = fastafile.split('.')[0]
    fastapath = 'gene_tree_files/fastas/{}'.format(fastafile)
    alnfile = './{}.aln'.format(fastabase)
    align = '{} {} > {}'.format(mafftbase,fastapath,alnfile)
    os.system(align)
    AlignIO.convert(alnfile,'fasta','gene_tree_files/nexus/{}.nex'.format(fastabase),'nexus',alphabet=Alphabet.generic_dna)
    os.remove(alnfile)
    return None

def parallelAlignFastas(processes):
    pool = ThreadPool(processes)
    fastas = os.listdir('gene_tree_files/fastas')
    alnfnc = partial(alignGeneFamilies)
    alnfiles = list(tqdm(pool.imap(alnfnc,fastas),total=len(fastas),desc='aligning'))
    return None

def buildTree(nexfile):
    nexbase = nexfile.split('.')[0]
    mbf = """set autoclose=yes nowarn=yes
execute gene_tree_files/nexus/{}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file={}
sump burnin=250
sumt burnin=250
quit""".format(nexfile,nexbase)
    mbscriptname = '{}_mrbayes_script.txt'.format(nexbase)
    open(mbscriptname,'w').write(mbf)
    logfilename = '{}_mrbayes_species_tree_log.txt'.format(nexbase)
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    os.system("/home/sid/thesis_SidReed/scripts/renameGeneTreeFiles.sh {}".format(nexbase))
#    treedir = 'gene_tree_files/trees/{}'.format(nexbase)
#    os.system('mkdir -p {}'.format(treedir))
#    for fp in glob.glob('./{}*'.format(nexbase)):
#        os.rename(fp,os.path.join(os.getcwd(),treedir,fp))
    return None

def parallelBuildTrees(processes):
    pool = ThreadPool(processes)
    alnfiles = os.listdir('gene_tree_files/nexus')
    treefnc = partial(buildTree)
    treefiles = list(tqdm(pool.imap(treefnc,alnfiles),total=len(alnfiles),desc='trees'))
    return None

def main(familys,nucFasta,minsize,processes):
    #set up out directories
    outdir = 'gene_tree_files'
    os.system('mkdir -p {}'.format(outdir))
    os.system('mkdir -p {}/fastas'.format(outdir))
    os.system('mkdir -p {}/nexus'.format(outdir))
    os.system('mkdir -p {}/trees'.format(outdir))
    #pick families
    family_idxs = pickGeneFamiliesForTrees(familys,minsize)
    print('Using {} of {} total families'.format(len(family_idxs),len(familys)))
    headerlists = {famidx:familys[famidx] for famidx in family_idxs}
    headerTuples = [(famidx,familys[famidx]) for famidx in family_idxs]
    #write families to fasta files
    seqs_to_write = extractSeqsParallel(headerTuples,nucFasta)
    seqs_to_write = {name:fixNames(seqs) for (name,seqs) in seqs_to_write}
    writeFastas(seqs_to_write)
    #align gene family fastas and convert to nexus
    parallelAlignFastas(processes)
    #build gene trees with mrbayes from alingments
    parallelBuildTrees(processes)
    return None


if __name__ == '__main__':
    famidx = json.load(open(sys.argv[1]))
    nucFasta = sys.argv[2]
    minsize = 0.4
    processes = 16
    main(famidx,nucFasta,minsize,processes)

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
