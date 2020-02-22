import os
import re
import sys
import json
import glob
from tqdm import tqdm
from Bio.Alphabet import IUPAC
from Bio import SeqIO,Seq,SeqRecord

#terms taht indicate a genes is annotated as a Mobile Genetic Element (MGE)
#any genes matching these terms (headers) are removed from the fastas
mobileElements = ['phage',
                  'mobile.element',
                  'insertion.sequence*',
                  'viral.element',
                  'transposase',
                  '[^sA-Z]IS\d',
                  '^IS\d',
                  'viral.([^A-Z]|enhancin)',
                  'holin']

def getHeaderAsDict(seqrecord,meregex):
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
    proteinname = [x for x in features if 'protein' in x][0].split('=')[1]
    notadd = any(rex.search(proteinname) for rex in meregex)
    for feat in features:
        if notadd:
            return None
        else:
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

def fixFasta(fasta,meregex,outname=-1):
    if outname == -1:
        outname = fasta
    allinfodirs = []
    fixedrecords = []
    for seqrecord in SeqIO.parse(open(fasta),'fasta'):
        infodir = getHeaderAsDict(seqrecord,meregex)
        if infodir != None:
            allinfodirs.append(infodir)
            fixedrecord = fixHeader(infodir)
            fixedrecords.append(fixedrecord)
        else:
            pass
    SeqIO.write(fixedrecords,outname,'fasta')
    return allinfodirs

def main(fasta_dir):
    allinfodirs = []
    meregex= [re.compile(x) for x in mobileElements]
    globstr = '{}/*.f*'.format(fasta_dir)
    for fasta in tqdm(glob.glob(globstr),desc='fixingfastas'):
        fastapath = os.path.join(os.getcwd(),fasta)
        infodir = fixFasta(fastapath,meregex)
        allinfodirs.extend(infodir)
    with open('fasta_headers_info.json','w') as outf:
        json.dump(allinfodirs,outf)
    return None


if __name__ == '__main__':
    fasta_dir  = sys.argv[1] #dir with fasta files /data/genus_data/$GENUS/nucleotide
    main(fasta_dir)
