import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

#ARGS
#arg 1 is the json with the list of gene families (including singletons verified by blasting against NR database)
#arg 2 should be path to a directory that contains all .faa files for every taxa in the matrix
## the gene names in the files in arg 2 should match the gene names from the families in arg 1
#arg 3 is the genus name
#python scripst/createPAMatrix.py gene_families.json protein_fastas/ genusname

def getOrganismGeneList(fastapath):
    genelist = []
    with open(fastapath) as fasta:
        for record in SeqIO.parse(fasta,'fasta'):
            genelist.append(record.id)
    return genelist

def allOrganismsGeneLists(fastadirpath):
    organismGenes = {}
    for fastapath in os.listdir(fastadirpath):
        if fastapath.endswith('.faa'):
            fullfastapath = os.path.join(fastadirpath,fastapath)
            organism = fastapath.split('.faa')[0]
            organismGenes[organism] = getOrganismGeneList(fullfastapath)
    return organismGenes

def buildmatrix(organismGenes,geneFamilies):
    pamat = np.zeros((len(list(organismGenes.keys())),len(list(geneFamilies.keys()))))
    orgIdxs = {i:o for i,o in enumerate(organismGenes.keys())}
    famIdxs = {i:f for i,f in enumerate(geneFamilies.keys())}
    for oidx,org in tqdm(orgIdxs.items(),total=len(list(orgIdxs.keys()))):
        for fidx,fam in tqdm(famIdxs.items(),total=len(list(famIdxs.keys()))):
            orgGenes = set(organismGenes[org])
            famGenes = set(geneFamilies[fam])
            pamat[oidx][fidx] = len(orgGenes.intersection(famGenes))
    return pamat,orgIdxs,famIdxs

def matrix_to_binary(pamat):
    binarize = lambda x:np.where(x > 0, 1, 0)
    vecbin = np.vectorize(binarize)
    return vecbin(pamat)

def namegenerator(basename):
#    matname = 'pa_matrix_{}.npy'.format(basename)
#    biname = 'binary_pa_matrix_{}.npy'.format(basename)
#    orgidxname = 'row_organism_idxs_{}.json'.format(basename)
#    famidxname = 'column_indexes_families_{}.json'.format(basename)
    matname = 'pa_matrix.npy'
    biname = 'binary_pa_matrix.npy'
    orgidxname = 'row_organism_idxs.json'
    famidxname = 'column_indexes_families.json'
    return matname,biname,orgidxname,famidxname

def main(families,fastadir):
    families = json.load(open(families))
    orggenes = allOrganismsGeneLists(fastadir)
    pamat = buildmatrix(orggenes,families)
    return pamat


if __name__ == '__main__':
    pamat,orgidxs,famidxs = main(sys.argv[1],sys.argv[2])
    matn,binn,orgidxn,famidxn = namegenerator(sys.argv[3])

    np.save(open(matn,'wb'),pamat)
    print('\nP/A matrix saved to {}'.format(matn))

    json.dump(orgidxs,open(orgidxn,'w'))
    print('row index organisms saved to {}'.format(orgidxn))
    json.dump(famidxs,open(famidxn,'w'))
    print('col index families saved to {}'.format(famidxn))
    print('family numbers in {} correspond to those in {}'.format(famidxn,sys.argv[1]))

    binmat = matrix_to_binary(pamat)
    np.save(open(binn,'wb'),binmat)
    print('binary P/A matrix saved to {}'.format(binn))

