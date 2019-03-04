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
#arg 1 is actualy genefamilies dict
#arg 2 is concatenated fasta fle with all genes (nuc seqs)
#arg 3 is genusname

def getNumTaxa():
    return len(glob.glob('./protein_fastas/*.faa'))

def toInclude(members,fam,maxtax,lower,upper):
    setsize = len(set([x.split(':')[0] for x in members]))
    setsizebool = (lower <= setsize) and (setsize <= upper)
    sizebool = (lower <= len(members)) and (len(members) <= upper)
    return (setsizebool and sizebool)

def pickGeneFamilies(genefamilies,maxtax,inclusion):
    #list of families that have exactly maxtax members and contain genes from every taxa being considered
    #return [fam for fam,members in tqdm(genefamilies.items(),total=len(genefamilies),desc='picking') if (len(set([x.split(':')[0] for x in members])) == maxtax) and (len(members) == maxtax)]
    #inclusion function
    if inclusion == 0:
        lower,upper = maxtaxx,maxtax
    elif (type(inclusion) == float) and (inclusion >= 0 and inclusion < 1):
        lower,upper = int(maxtax*(1-inclusion)),int(maxtax*(1+inclusion))
    elif (type(inclusion) == int):
        lower,upper = maxtax-inclusion,maxtax+inclusion
    else:
        raise ValueError('invalid inclusion criteria')
    print('trees for species tree must have between {} and {} leaves (total tax is {})'.format(lower,upper,maxtax))
    inclusionfnc = partial(toInclude,maxtax=maxtax,lower=lower,upper=upper)
    return [fam for fam,members in tqdm(genefamilies.items(),total=len(genefamilies),desc='picking') if inclusionfnc(members,fam)]

#def iterativePickFamilies(families,maxtax,inclusion,iterations=5,threshold=50):
#    enoughFamilies = False
#    while not enoughFamilies:
#        famlist = pickGeneFamilies(families,maxtax,inclusion)
#        if len(famlist)

def extractSeqsForTree(headerlist,nucFasta,outdir,genefam):
    seqRecordlist = []
    seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
    seqRecordlist = [seq for seq in seqiter if seq.id in headerlist] #put sequences in list
    fastafilename = '{}/fastas/{}.fna'.format(outdir,genefam)
    return fastafilename, seqRecordlist

def fixNames(seqlist):
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

def concatNexAlns(outdir,genusname):
    """Combine multiple nexus data matrices in one partitioned file.
    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this """
    nexdir = '{}/nexus/'.format(outdir)
    filelist = [x for x in os.listdir(nexdir) if x.endswith('.nex')]
    nexi = [(os.path.join(nexdir,fname), Nexus.Nexus(os.path.join(nexdir,fname))) for fname in filelist]
    coutname = '{}/concat_aln_species_tree_{}.nex'.format(outdir,genusname)
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open(coutname,'w'))
    return coutname

def buildSpeciesTree(concat_nexus_aln,outdir,genusname):
    mbf =\
"""set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file=tmp_species_tree_{}
sump burnin=250
sumt burnin=250
quit""".format(concat_nexus_aln,genusname)
    mbscriptname = 'tmp_mrbayes_script.txt'
    open(mbscriptname,'w').write(mbf)
    logfilename = 'tmp_mrbayes_species_tree_log.txt'
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    treedir= '{}/species_tree_{}'.format(outdir,genusname)
    os.system('mkdir -p {}'.format(treedir))
    os.rename(mbscriptname,os.path.join(os.getcwd(),treedir,mbscriptname))
    for fp in glob.glob('./*tmp_*'):
        os.rename(fp,os.path.join(os.getcwd(),treedir,fp.replace('tmp_','')))
    return treedir

def main(familys,nucFasta,genusname,genes,processes,inclusion):
    #set up out directories
    outdir = 'species_tree_files'
    os.system('mkdir -p {}'.format(outdir))
    os.system('mkdir -p {}/fastas'.format(outdir))
    os.system('mkdir -p {}/nexus'.format(outdir))
    print('picking gene families for species tree...')
    maxtax = getNumTaxa()
    family_idxs = pickGeneFamilies(familys,maxtax,inclusion)
    print('Using {} of {} qualified families of {} total families'.format(genes,len(family_idxs),len(familys)))
    headerlists = {famidx:familys[famidx] for famidx in family_idxs[:genes]}
    print('writing gene families to fastas...')
    seqs_to_write = [extractSeqsForTree(heads,nucFasta,outdir,fam) for fam,heads in tqdm(headerlists.items(),total=len(list(headerlists.keys())),desc='getseqs')]
    seqs_to_write = {name:fixNames(seqs) for (name,seqs) in seqs_to_write}
    writeFastas(seqs_to_write)
    print('aligning and concatenating genes...')
    parallelAlignFastas(outdir,processes)
    concatfilename = concatNexAlns(outdir,genusname)
    print('building species tree with concatenated alingment...')
    treedir = buildSpeciesTree(concatfilename,outdir,genusname)
    return None


if __name__ == '__main__':
    genefamilies = json.load(open(sys.argv[1],'r'))
    nucFasta = sys.argv[2]
    genusname = sys.argv[3]
    genes_for_species_tree = 50
    processes = 32
    inclusion = 2
    main(genefamilies,
         nucFasta,
         genusname,
         genes_for_species_tree,
         processes,
         inclusion)

#DEPRECIATED
def matrix_to_binary(pamat):
    binarize = lambda x:np.where(x > 0, 1, 0)
    vecbin = np.vectorize(binarize)
    return vecbin(pamat)

def pickGeneFamilies(pamat,columnindex,genes):
    bmat = matrix_to_binary(pamat)
    sizes = np.apply_along_axis(sum,0,bmat)
    #family_idxs = [str(i) for i,pavec in enumerate(pamat.T) if np.array_equal(pavec,np.ones(pamat.shape[0]))]
    family_idxs = [str(i) for i,sumvec in enumerate(sizes) if sumvec == bmat.shape[0]]
    print('using {} of {} usable of {} total genes for species tree'.format(genes,len(family_idxs),pamat.shape[1]))
    return family_idxs[:genes]

