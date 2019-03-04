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
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO,AlignIO,Phylo,Alphabet
from multiprocessing.dummy import Pool as ThreadPool

#will extract all genes present in only 1 copy in every taxa and spit to a fasta file for alignment and tree creation
#args
#arg 1 is concatenated fasta fle with all genes (nuc seqs)
#arg 2 is fasta info json

#GET ALL FAMILIES THAT CONTAIN 16S, Align, combine, build tree

def get16SHeaders(fastaInfoJson):
    return [nucDir['header'] for nucDir in tqdm(fastaInfoJson,desc='picking') if '16S rRNA' in nucDir['protein']]

def get16SFamilies(headers16, geneFamilies):
    fams16 = {}
    print(len(geneFamilies))
    for fam,genes in geneFamilies.items():
        if len(set(genes).intersection(set(headers16))) != 0:
            fams16[fam] = genes
    return fams16

def fixNames(seqlist):#list of seqrecords
    fixedseqs = []
    found_taxa = []
    for cseq in seqlist:
        new_name = cseq.name.split(':')[0]
        if not new_name in found_taxa:
            fixedseqs.append(SeqRecord(cseq.seq,
                                       name=new_name,
                                       id=new_name,
                                       description=''))
            found_taxa.append(new_name)
    return fixedseqs

def extractSeqsForTree(famdict,nucFasta):
    namesAndRecords = []
    for fam,headerlist in famdict.items():
        seqRecordlist = []
        seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
        seqRecordlist = [seq for seq in seqiter if seq.id in headerlist]
        seqRecordlist = fixNames(seqRecordlist)
        fastafilename = 'species_tree_files/fastas/{}.fna'.format(fam)
        namesAndRecords.append((fastafilename,seqRecordlist))
    return namesAndRecords

def writeFastas(nameseqiter):
    for (name,seqs) in nameseqiter:
        with open(name,'w') as outfile:
            SeqIO.write(seqs,outfile,'fasta')#write all seqs to a fasta file
    return None

def alignGeneFamilies(fastafile):
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    fastabase = fastafile.split('.')[0]
    fastaname = 'species_tree_files/fastas/{}'.format(fastafile)
    alnfile = 'species_tree_files/tmp_{}.aln'.format(fastabase)
    align = '{} {} > {}'.format(mafftbase,fastaname,alnfile)
    os.system(align)
    nexfile = 'species_tree_files/nexus/{}.nex'.format(fastabase)
    AlignIO.convert(alnfile,'fasta',nexfile,'nexus',alphabet=Alphabet.generic_dna)
    os.remove(alnfile)
    return None

def parallelAlignFastas(processes):
    pool = ThreadPool(processes)
    fastas = os.listdir('species_tree_files/fastas')
    alnfiles = list(tqdm(pool.imap(alignGeneFamilies,fastas),total=len(fastas),desc='aligning'))
    return None

def concatNexAlns():
    """Combine multiple nexus data matrices in one partitioned file.
    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this """
    nexdir = 'species_tree_files/nexus/'
    filelist = [x for x in os.listdir(nexdir) if x.endswith('.nex')]
    nexi = [(os.path.join(nexdir,fname), Nexus.Nexus(os.path.join(nexdir,fname))) for fname in filelist]
    coutname = 'species_tree_files/concat_aln_species_tree.nex'
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open(coutname,'w'))
    return coutname

def buildSpeciesTree(concat_nexus_aln):
    mbf =\
"""set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file=tmp_species_tree
sump burnin=250
sumt burnin=250
quit""".format(concat_nexus_aln)
    mbscriptname = 'tmp_mrbayes_script.txt'
    open(mbscriptname,'w').write(mbf)
    logfilename = 'tmp_mrbayes_species_tree_log.txt'
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    treedir= 'species_tree_files/species_tree'
    os.system('mkdir -p {}'.format(treedir))
    os.rename(mbscriptname,os.path.join(os.getcwd(),treedir,mbscriptname))
    for fp in glob.glob('./*tmp_*'):
        os.rename(fp,os.path.join(os.getcwd(),treedir,fp.replace('tmp_','')))
    return treedir

def main(allNucFasta,fastaInfoJson,geneFamilies,processes,maxfams):
    base = 'species_tree_files'
    os.system('mkdir -p {}'.format(base))
    os.system('mkdir -p {}/fastas/'.format(base))
    os.system('mkdir -p {}/nexus'.format(base))
    headers16s = get16SHeaders(fastaInfoJson)
    print(len(headers16s))
    fams16 = get16SFamilies(headers16s,geneFamilies)
    print(len(fams16))
    namesAndRecords = extractSeqsForTree(fams16,allNucFasta)[:maxfams]
    print(len(namesAndRecords))
    if maxfams > len(fams16):
        print('using {} of {} available 16S rRNA fams of {} total fams for species tree'.format(len(fams16),len(fams16),len(list(geneFamilies.keys()))))
    else:
        print('using {} of {} available 16S rRNA fams of {} total fams for species tree'.format(maxfams,len(fams16),len(list(geneFamilies.keys()))))
    writeFastas(namesAndRecords)
    parallelAlignFastas(processes)
    combinedNexusName = concatNexAlns()
    buildSpeciesTree(combinedNexusName)
    return None


if __name__ == '__main__':
    nucFasta = sys.argv[1]
    fastaInfoJson = json.load(open(sys.argv[2]))
    geneFamilies = json.load(open(sys.argv[3],'r'))
    maxfams = 50
    processes = 32
    main(nucFasta,
         fastaInfoJson,
         geneFamilies,
         processes,
         maxfams,
         )

