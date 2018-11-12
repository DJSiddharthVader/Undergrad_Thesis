import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

#arg 1 is PA matrix (not binary)
#arg 2 is row oranisms index dict
#arg 3 is col families index dict
#arg 4 is actualy genefamilies dict
#arg 5 is concatenated fasta fle with all genes
#arg 6 is outname, no extension

#will extract all genes present in only 1 copy in every taxa and spit to a fasta file for alignment and tree creation

def getGeneFamilies(pamat,columnindex):
    familylist = []
    for family in pamat.T:
        if sum(family) == np.ones(len(family)):
            familylist.append(family)
    return familylist

def familyToList(familylist,familyindex,fams=50):
#list  of lists, each is a genefamily with 1 member in every organism
    return [familyindex[colidx] for colidx in familylist]

def listToFasta(genelist,fastaname,outname):
    singlerecords = []
    with open(fastaname) as fasta:
        for record in SeqIO.parse(fasta,'fasta'):
            if record.id in genelist:
                singlerecords.append(record)
    foutname = 'singletons_{}.faa'.format(outname)
    with open(foutname,'w') as singleout:
        SeqIO.write(singlerecords,foutname,'fasta')
    return foutname


if __name__ == '__main__':
    print(None)
