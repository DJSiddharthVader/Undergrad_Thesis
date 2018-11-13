import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Nexus import Nexus
from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

#arg 1 is PA matrix (not binary)
#arg 2 is col families index dict
#arg 3 is actualy genefamilies dict
#arg 4 is concatenated fasta fle with all genes (nuc seqs)
#arg 5 is outname, no extension

#will extract all genes present in only 1 copy in every taxa and spit to a fasta file for alignment and tree creation

def getGeneFamilies(pamat,columnindex):
    familyidxlist = []
    print(pamat.T.shape)
    orgs = pamat.T.shape[1]
    for i,family in enumerate(pamat.T):
        if np.array_equal(family,np.ones(orgs)):
            familyidxlist.append(columnindex[str(i)])
    return familyidxlist #gene family colidxs with only 1 gene in every taxa

def familyToList(familylist,familyindex,fams):
    #list  of lists, each is a genefamily with 1 member in every organism
    return [(colidx,familyindex[colidx]) for colidx in familylist][:fams] #elements are tuples, fam index and gene list

def listToFasta(genelist,fastaname,outdir):
    (name,genelist) = genelist
    singlerecords = []
    with open(fastaname) as fasta:
        for record in SeqIO.parse(fasta,'fasta'):
            if record.id in genelist:
                singlerecords.append(record)
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
        return outdir
    pool = ThreadPool(processes)
    parlistToFasta = partial(listToFasta,fastaname=nucFasta,outdir=outdir)
    for _ in tqdm(pool.imap_unordered(parlistToFasta,familylist),total=len(familylist),desc='extracting families'):
        pass
    return outdir

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
    print(alltaxa)
    print(len(alltaxa))
    print([len(x[1].taxlabels) for x in nexi])
    missingTaxa = [False if set(x[1].taxlabels) == alltaxa else True for x in nexi]
    print(set(missingTaxa))
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

def main(pamat,columnindex,familyindex,nucFasta,outname,processes,fams):
    fastadir = 'species_tree_fastas_{}'.format(outname)
    fastadir = writeAllSpeciesTreeFastas(pamat,columnindex,familyindex,nucFasta,fastadir,processes,fams)
    print('fasta files for species tree genes written to {}'.format(fastadir))
    nexDir = nexusAlnSTreeFastas(fastadir,fams)
    print('fasta alingments for species tree written to {}'.format(fastadir))
    print('nexus alingments for species tree written to {}'.format(nexDir))
    coutname = concatNexAlns(nexDir,outname,same_taxa=True)
    print('concat nexus aln file for mrbayes written to {}'.format(coutname))

if __name__ == '__main__':
    pamat = np.load(sys.argv[1])
    columnindex = json.load(open(sys.argv[2]))
    familyindex = json.load(open(sys.argv[3]))
    nucFasta = sys.argv[4]
    outname = sys.argv[5]
    processes = 8
    fams = 50
    main(pamat,columnindex,familyindex,nucFasta,outname,processes,fams)


def checkTaxa(nexi):#from https://biopython.org/wiki/Concatenate_nexus
    """Verify Nexus instances have the same taxa information.
    Checks that nexus instances in a list [(name, instance)...] have
    the same taxa, provides useful error if not and returns None if
    everything matches"""
    first_taxa = nexi[0][1].taxlabels
    for name, matrix in nexi[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]
        if first_only:
            missing = ', '.join(first_only)
            msg='{} taxa {} not in martix {}'.format(nexi[0][0],missing,name)
            raise Nexus.NexusError(msg)
        elif new_only:
            missing = ', '.join(new_only)
            msg = '{} taxa {} not in all matrices'.format(name, missing)
            raise Nexus.NexusError(msg)
    return None # will only get here if it hasn't thrown an exception
