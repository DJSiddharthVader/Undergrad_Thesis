import os
import sys
import json
import glob
import string
from tqdm import tqdm
from Bio.Alphabet import IUPAC
from Bio import SeqIO,Seq,SeqRecord

def getHeaderAsDict(seqrecord):
    raw = seqrecord.id.strip('lcl|')
    acc = '_'.join(raw.split('_')[:2])
    wp  = '_'.join(raw.split('_')[3:5])
    infodir = {'organism_accession':acc,
               'protein_accession':wp,
               'sequence':str(seqrecord.seq),
               'header':'{}:{}'.format(acc,wp),
               }
    desc = ' '.join(seqrecord.description.split(' ')[1:]).strip('[]')
    features = desc.split('] [')
    for feat in features:
        infodir[feat.split('=')[0]] = feat.split('=')[1]
    return infodir

def strToSeqObj(seq):
    aas = set(IUPAC.protein.letters)
    if len(set(seq).intersection(aas)) != 0:
        seqobj = Seq.Seq(seq,IUPAC.protein)
    else:
        seqobj = Seq.Seq(seq,IUPAC.unambiguous_dna)
    return seqobj

def fixHeader(infodir):
    newid = "{}:{}".format(infodir['organism_accession'],infodir['protein_accession'])
    seqobj = strToSeqObj(infodir['sequence'])
    fixedrecord = SeqRecord.SeqRecord(seqobj,
                                    id=newid,
                                    name=newid,
                                    description='')
    return fixedrecord

def fixFasta(fasta,outname=-1):
    if outname == -1:
        outname = fasta
    allinfodirs = []
    fixedrecords = []
    for seqrecord in SeqIO.parse(open(fasta),'fasta'):
        infodir = getHeaderAsDict(seqrecord)
        allinfodirs.append(infodir)
        fixedrecord = fixHeader(infodir)
        fixedrecords.append(fixedrecord)
    SeqIO.write(fixedrecords,outname,'fasta')
    return allinfodirs

def main(fasta_dir,genus):
    allinfodirs = []
    globstr = '{}/*.f*a'.format(fasta_dir)
    for fasta in tqdm(glob.glob(globstr),desc='fixingfastas'):
        fastapath = os.path.join(os.getcwd(),fasta)
        infodir = fixFasta(fastapath)
        allinfodirs.extend(infodir)
    with open('{}_fastas_info.json'.format(genus),'w') as outf:
        json.dump(allinfodirs,outf)
    return None

if __name__ == '__main__':
    fasta_dir  = sys.argv[1]
    genus_name = sys.argv[2]
    main(fasta_dir,genus_name)
