import os
import sys
import json
import subprocess
from tqdm import tqdm
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO,AlignIO,Alphabet

def get16SHeaders(fastaInfoJson):
    return [nucDir['header'] for nucDir in tqdm(fastaInfoJson,desc='picking') if '16S rRNA' in nucDir['protein']]

def fixRecord(seqlist):#list of seqrecords
    fixedseqs = []
    found_taxa = []
    for cseq in seqlist:
        new_name = cseq.name.split(':')[0]
        if not new_name in found_taxa:
            found_taxa.append(new_name)
            fixedseqs.append(SeqRecord(cseq.seq,
                                       name=new_name,
                                       id=new_name,
                                       description=''))
    return fixedseqs

def extractSeqsForTree(headerlist,nucFasta):
    seqRecordlist = []
    seqiter = list(SeqIO.parse(open(nucFasta),'fasta'))
    seqRecordlist = [seq for seq in seqiter if seq.name in headerlist]
    seqRecordlist = fixRecord(seqRecordlist)
    fastafile = 'species_tree_files/16s_rRNA_seqs.fna'
    SeqIO.write(seqRecordlist,fastafile,'fasta')
    return fastafile

def alignFasta(fastafile):
    dirname,filename = fastafile.split('/')
    mafftbase = 'mafft --quiet --localpair --maxiterate 1000'
    outfile = '{}/{}.aln'.format(dirname,filename.split('.')[0])
    cmd = '{} {} > {}'.format(mafftbase,fastafile,outfile)
    os.system(cmd)
    nexname = outfile.split('.')[0] + '.nex'
    AlignIO.convert(outfile,'fasta',nexname,'nexus',alphabet=Alphabet.IUPAC.unambiguous_dna)
    return nexname

def buildSpeciesTree(nexus_aln):
    mbf =\
"""set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file=species_tree_files/species_tree
sump burnin=250
sumt burnin=250
quit""".format(nexus_aln)
    mbscriptname = 'species_tree_files/mrbayes_script.txt'
    open(mbscriptname,'w').write(mbf)
    logfilename = 'species_tree_files/mrbayes_species_tree_log.txt'
    with open(mbscriptname,'rb',0) as mbscript, open(logfilename,'wb',0) as logfile:
        logtxt = subprocess.run(['mb'],
                                stdin=mbscript,
                                stdout=logfile,
                                check=True)
    return None

def main(fasta,info,genefamilies,maxfams):
    os.system('mkdir -p species_tree_files')
    headers16s = get16SHeaders(info)
    if maxfams > len(headers16s):
        print('using {} of {} available 16S rRNA genes for species tree'.format(len(headers16s),len(headers16s)))
    else:
        print('using {} of {} available 16S rRNA genes for species tree'.format(maxfams,len(headers16s)))
    print('getting sequences')
    fastafile = extractSeqsForTree(headers16s,fasta)
    print('aligning sequences')
    alnfile = alignFasta(fastafile)
    print('building species tree')
    buildSpeciesTree(alnfile)
    return None


if __name__ == '__main__':
    nucFasta = sys.argv[1]
    fastaInfoJson = json.load(open(sys.argv[2]))
    geneFamilies = json.load(open(sys.argv[3],'r'))
    maxfams = 50
    main(nucFasta,
         fastaInfoJson,
         geneFamilies,
         maxfams,
         )
