import sys
from tqdm import tqdm
from Bio import SeqIO

def dedupFasta(fasta):
    uniquerecords = []
    uniqueseqs = []
    for seqrecord in tqdm(SeqIO.parse(open(fasta),'fasta')):
        if str(seqrecord.seq) not in uniqueseqs:
            uniquerecords.append(seqrecord)
            uniqueseqs.append(str(seqrecord.seq))
    return uniquerecords

def main(fasta,outname=-1):
    if outname == -1:
        outname = fasta
    deduped = dedupFasta(fasta)
    SeqIO.write(deduped,fasta,'fasta')
    return None

if __name__ == '__main__':
    main(sys.argv[1])

