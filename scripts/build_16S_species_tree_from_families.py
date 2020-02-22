import os
import re
import sys
import json
import glob
import subprocess
from tqdm import tqdm
from Bio.Nexus import Nexus
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO,Alphabet,SeqIO
from multiprocessing.dummy import Pool as ThreadPool

#will extract all genes present in only 1 copy in every taxa and spit to a fasta file for alignment and tree creation
#args
#arg 1 is concatenated fasta fle with all genes (nuc seqs)
#arg 2 is fasta info json

#GET ALL FAMILIES THAT CONTAIN 16S, Align, combine, build tree
base_dir = 'species_tree_16S'

def getNumTaxa():
    return len(glob.glob('./protein_fastas/*.faa'))

def get16SHeaders(fastaInfoJson):
    return [nucDir for nucDir in tqdm(fastaInfoJson,desc='picking') if '16S rRNA' in nucDir['protein']]

def pick16Sprots(headers16s):
    headerlists = []
    numtaxa = getNumTaxa()
    for p in set([x['protein'] for x in headers16s]):
        ponly = [x for x in headers16s if x['protein'] == p]
        if len(set([x['organism_accession'] for x in ponly])) == numtaxa:
            headerlists.append(ponly)
    return headerlists

def getFamDict16s(headerlists):
    famdict = {}
    for hlist in headerlists:
        famdict[hlist[0]['protein']] = [x['header'] for x in hlist]
    return famdict

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
    for name,headerlist in famdict.items():
        seqRecordlist = []
        seqiter = SeqIO.parse(open(nucFasta),'fasta') #open fasta file
        seqRecordlist = [seq for seq in seqiter if seq.id in headerlist]
        seqRecordlist = fixNames(seqRecordlist)
        name = re.sub('[^0-9a-zA-Z]+','-',name)
        fastafilename = '{}/fastas/{}.fna'.format(base_dir,name)
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
    fastaname = '{}/fastas/{}'.format(base_dir,fastafile)
    alnfile = 'tmp_{}.aln'.format(fastabase)
    align = '{} {} > {}'.format(mafftbase,fastaname,alnfile)
    os.system(align)
    nexfile = '{}/nexus/{}.nex'.format(base_dir,fastabase)
    AlignIO.convert(alnfile,'fasta',nexfile,'nexus',alphabet=Alphabet.generic_dna)
    os.remove(alnfile)
    return None

def parallelAlignFastas(processes):
    pool = ThreadPool(processes)
    fastas = os.listdir('{}/fastas'.format(base_dir))
    alnfiles = list(tqdm(pool.imap(alignGeneFamilies,fastas),total=len(fastas),desc='aligning'))
    return None

def concatNexAlns():
    """Combine multiple nexus data matrices in one partitioned file.
    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this """
    nexdir = '{}/nexus/'.format(base_dir)
    filelist = [x for x in os.listdir(nexdir) if x.endswith('.nex')]
    nexi = [(os.path.join(nexdir,fname), Nexus.Nexus(os.path.join(nexdir,fname))) for fname in filelist]
    coutname = '{}/concat_aln_species_tree.nex'.format(base_dir)
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
    treedir= '{}/species_tree'.format(base_dir)
    os.system('mkdir -p {}'.format(treedir))
    os.rename(mbscriptname,os.path.join(os.getcwd(),treedir,mbscriptname))
    for fp in glob.glob('./*tmp_*'):
        os.rename(fp,os.path.join(os.getcwd(),treedir,fp.replace('tmp_','')))
    return treedir


def main(allNucFasta,fastaInfoJson,processes,maxfams):
    #make dirs
    os.system('mkdir -p {}'.format(base_dir))
    os.system('mkdir -p {}/fastas/'.format(base_dir))
    os.system('mkdir -p {}/nexus'.format(base_dir))
    #Get 16S Gene Fastas
    headers16s = get16SHeaders(fastaInfoJson)
    fams16 = pick16Sprots(headers16s)
    famsdict16 = getFamDict16s(fams16[:maxfams])
    namesAndRecords = extractSeqsForTree(famsdict16,allNucFasta)
    writeFastas(namesAndRecords)
    #print how many fams being used
    if len(famsdict16) == 0:
        print('No 16S sequences common to all organisms')
        sys.exit(1)
    if maxfams > len(fams16):
        print('using {} of {} available 16S rRNA genes for species tree'.format(len(list(famsdict16.keys())),len(fams16)))
    else:
        print('using {} of {} available 16S rRNA genes for species tree'.format(maxfams,len(fams16)))
    print('aligning seqs')
    parallelAlignFastas(processes)
    combinedNexusName = concatNexAlns()
    print('running mrbayes')
    buildSpeciesTree(combinedNexusName)
    return None


if __name__ == '__main__':
    nucFasta = sys.argv[1]
    fastaInfoJson = json.load(open(sys.argv[2]))
    maxfams = 50
    processes = 16
    main(nucFasta,
         fastaInfoJson,
         processes,
         maxfams)

