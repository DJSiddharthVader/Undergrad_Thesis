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

def getGeneFamilies(pamat,columnindex):
    orgs = pamat.T.shape[1]
    return [columnindex[str(i)] for i,family in enumerate(pamat.T) if np.array_equal(family,np.ones(orgs))]
#    familyidxlist = []
#    for i,family in enumerate(pamat.T):
#        if np.array_equal(family,np.ones(orgs)):
#            familyidxlist.append(columnindex[str(i)])
#    return familyidxlist #gene family colidxs with only 1 gene in every taxa

def familyToList(familylist,familyindex,fams):
    #list  of lists, each is a genefamily with 1 member in every organism
    #elements are tuples, fam index and gene list
    return [(colidx,familyindex[colidx]) for colidx in familylist][:fams]

def listToFasta(genelist,fastaname,outdir):
    (name,genelist) = genelist
    singlerecords = []
    with open(fastaname) as fasta:
        for record in SeqIO.parse(fasta,'fasta'):
            if (record.id in genelist) and (str(record.seq) != ''):
                newrecord = record
                newrecord.id = record.name.split(':')[0].strip('>')
                newrecord.name = record.name.split(':')[0].strip('>')
                singlerecords.append(newrecord)
    foutname = 'gene_{}.fna'.format(name)
    absoutname = os.path.join(outdir,foutname)
    with open(absoutname,'w') as singleout:
        SeqIO.write(singlerecords,absoutname,'fasta')
    return absoutname

def writeAllSpeciesTreeFastas(pamat,columnindex,familyindex,nucFasta,outdir,processes,fams):
    os.system('mkdir -p {}'.format(outdir))
    familyidxlist = getGeneFamilies(pamat,columnindex)
    familylist = familyToList(familyidxlist,familyindex,fams)
    if len(os.listdir(outdir)) == fams: #previously cached results
        #return outdir
        pass
    pool = ThreadPool(processes)
    parlistToFasta = partial(listToFasta,fastaname=nucFasta,outdir=outdir)
    for _ in tqdm(pool.imap_unordered(parlistToFasta,familylist),total=len(familylist),desc='extracting families'):
        pass
    return outdir,len(familylist)

def nexusAlnSTreeFastas(sTreeDir,fams):
    faadir = 'aligned_faa_{}'.format(sTreeDir)
    os.system('mkdir -p {}'.format(faadir))
    if len(os.listdir(faadir)) != fams:
        mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
        for fasta in tqdm(os.listdir(sTreeDir),desc='alining'):
            if fasta.endswith('.fna'):
                fout = fasta.strip('.fna')
                mafft ='{} {}/{} > {}/{}.aln'.format(mafftbase,sTreeDir,fasta,faadir,fout)
                os.system(mafft)
    nexDir = 'aligned_nex_{}'.format(sTreeDir)
    os.system('mkdir -p {}'.format(nexDir))
    if len(os.listdir(nexDir)) != fams:
        for aln in tqdm(os.listdir(faadir),desc='to nexus'):
            AlignIO.convert('{}/{}'.format(faadir,aln),'fasta','{}/{}.nex'.format(nexDir,aln.strip('.aln')),'nexus',alphabet=Alphabet.generic_dna)
    return nexDir

def checkTaxa(nexi):
    alltaxa = set([y for x in nexi for y in x[1].taxlabels])
    missingTaxa = [False if set(x[1].taxlabels) == alltaxa else True for x in nexi]
    return False

def concatNexAlns(nexDir,outname,same_taxa=True):#from https://biopython.org/wiki/Concatenate_nexus
    """Combine multiple nexus data matrices in one partitioned file.
    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this """
    filelist = [x for x in os.listdir(nexDir) if x.endswith('.nex')]
    nexi = [(os.path.join(nexDir,fname), Nexus.Nexus(os.path.join(nexDir,fname))) for fname in filelist]
    coutname = 'concat_stree_aln_{}.nex'.format(outname)
    if same_taxa:
        if not checkTaxa(nexi):
            combined = Nexus.combine(nexi)
            combined.write_nexus_data(filename=open(coutname,'w'))
            return coutname
    else:
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
    columnindex = json.load(open(sys.argv[2]))
    familyindex = json.load(open(sys.argv[3]))
    nucFasta = sys.argv[4]
    outname = sys.argv[5]
    processes = 8
    fams = 50
    main(pamat,columnindex,familyindex,nucFasta,outname,processes,fams)

