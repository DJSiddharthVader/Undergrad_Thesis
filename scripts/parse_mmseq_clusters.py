import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from functools import partial
from collections import defaultdict
from multiprocessing.dummy import Pool as ThreadPool

def chunklines(cfile,chunksize=100):
    alllines = open(cfile).readlines()
    numchunks = len(alllines)//chunksize + 1
    chunks = [alllines[i*chunksize:(i+1)*chunksize] for i in range(0,numchunks)]
    return chunks

def parse(lines):
    clusters = defaultdict(list)
    for line in lines:
        cluster,member = line.strip('\n').split('\t')
        clusters[cluster].append(member)
    return clusters

def mergeDOLs(dols):
    base = dols[0]
    singletons = []
    for d in tqdm(dols[1:]):
        for k,v in d.items():
            if len(v) == 1:
                singletons.extend(v)
            else:
                base[k].extend(v)
    finalcopy = defaultdict(list)
    for i, (k,v) in enumerate(base.items()):
        finalcopy['fam{}'.format(i)].extend(v)
    return finalcopy,singletons

def parseParallel(cfile,processes):
    pool = ThreadPool(processes)
    lines = chunklines(cfile)
    clusterlists = list(tqdm(pool.imap(parse,lines),total=len(lines),desc='parsing'))
    clusters, singletons = mergeDOLs(clusterlists)
    return clusters,singletons

def writeSingletons(genelist,allfasta):
    parser = SeqIO.parse(open(allfasta),'fasta')
    singlerecords = [seq for seq in parser if seq.id in genelist]
    SeqIO.write(singlerecords,'singletons.faa','fasta')
    return None

def main(clusters,allfasta,processes=32):
    clusters,singletons = parseParallel(clusters,processes)
    writeSingletons(singletons,allfasta)
    json.dump(clusters,open('gene_families.json','w'))

if __name__ == '__main__':
    #arg1 is clusters.tsv output, arg2 is all_proteins.faa file
    main(sys.argv[1],sys.argv[2])
