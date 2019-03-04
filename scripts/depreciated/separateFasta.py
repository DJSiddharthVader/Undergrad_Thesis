import os
import sys
from Bio import SeqIO

#splits fasta with nuc and prot into 2 fasta by type
def filtermixed(mixed):
    base = os.path.basename(mixed).split('.')[0]
    nucname = base + '.fna'
    protname = base + '.faa'
    nuclist,protlist = [],[]
    for record in SeqIO.parse(open(mixed),'fasta'):
        if all(c.upper() in 'ATCG' for c in record.seq): #DNA
            protlist.append(record)
        else:
            nuclist.append(record)
    with open(protname,'w') as pout:
        SeqIO.write(protlist,pout,'fasta')
    with open(nucname,'w') as nout:
        SeqIO.write(nuclist,nout,'fasta')
    return protname,nucname

if __name__ == '__main__':
    n1,n2 = filtermixed(sys.argv[1])
    print(n1)
    print(n2)
