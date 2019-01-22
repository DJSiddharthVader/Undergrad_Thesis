import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import networkx as nx
from Bio import Phylo
from ete3 import Tree
from tqdm import tqdm
import dendropy as dp
from functools import partial
from io import StringIO
from multiprocessing.dummy import Pool as ThreadPool

#Docs
#STEPS
#1. convert species nexus_tree to newick tree
#2. convert gene nexus_trees to newick trees
#3. split networks into subsets (how yet not sure)
#4. run HiDe on each subset to produce 1 network
#5. parse raw HiDe output to some computer readable format
#6. write to json file for use in networkx/igraph
#produces file direcories as follows
#network_files/
#    all_newick_trees/species.newick
#                     fam1.newick
#                     ...
#    subsets/
#        subset_01/species.newick
#                  fam1.newick
#                  fam2.newick
#                  fam4.newick
#                  ...
#        subset_02/species.newick
#                  fam3.newick
#                  fam2.newick
#                  fam6.newick
#                  ...
#        ...
#    networks_raw/all_gene_trees.txt
#                 subset_01.txt
#                 ...
#    networks_json/all_gene_trees.json
#                  subset_01.json
#                  ...

#make network
def writeNewickWithOutBrLen(treefile,outpath):
    nextree = dp.Tree.get(path=treefile,schema='nexus')
    treestr = nextree.as_string(schema='newick')
    treestr = re.sub(':\d\.[0-9\-e]+','',treestr)
    treestr = re.sub('\[\&[UR]\]','',treestr)
    if re.search('copy',treestr) != None:
        hascopy = 1
    else:
        hascopy = 0
        nwktree = Tree(treestr,format=0)
        nwktree.resolve_polytomy(recursive=True)
        #r = nwktree.get_midpoint_outgroup()
        #nwktree.set_outgroup(r)
        nwktree.write(format=9,outfile=outpath)
    return hascopy

def writeNewickParallel(treefile):
    basename = os.path.basename(treefile).split('.')[0]
    nwkfile = 'network_files/all_newick_trees/{}.newick'.format(basename)
    return writeNewickWithOutBrLen(treefile,nwkfile)

def fixNewickTrees(speciesdir,genedir,processes):
    #set up dirs
    maindir = 'network_files'
    os.system('mkdir -p {}/all_newick_trees'.format(maindir))
    #convert species tree to newick
    speciestreefile = glob.glob('{}/*.con.tre'.format(speciesdir))[0]
    writeNewickWithOutBrLen(speciestreefile,'{}/all_newick_trees/species.newick'.format(maindir))
    #convert gene trees to newick
    pool = ThreadPool(processes)
    treefiles = list(glob.iglob('{}/**/*.con.tre'.format(genedir)))
    hascopy = sum(list(tqdm(pool.imap(writeNewickParallel,treefiles),total=len(treefiles),desc='fixtrees')))
    print('{} of {} gene trees had duplicate taxa, not considered'.format(hascopy,len(os.listdir(genedir))))
    numgenetrees = len(os.listdir('network_files/all_newick_trees')) -1 #-1 for the species tree
    return '{}/all_newick_trees'.format(maindir), numgenetrees

#cd to hide dir
def runHideOnTreeDir(treedir,networkfile):
#run hide and capture output
#write network to file
#parse network
#write as graphml
    networktxt = subprocess.run(['./lua','score.lua',treedir],
                            stdout=networkfile,
                            check=True)
    return None

#parse network
def parseEdge(edge):
    membs = edge.split(' ')
    source = re.sub('\d*\-','',membs[0])
    sink = re.sub('\d*\-','',membs[2])
    isInternal = True
    if re.search('^LCA\(',source) == None:
        if re.search('^LCA\(',sink) == None:
            isInternal = False
    return source,sink,isInternal

def annotateCRISPR(network,crisprdata):
    networkcrisprdata = crisprdata[crisprdata['Accession Number'].isin(network.nodes())]
    crisprlabels = {row[1]['Accession Number']:\
                    {header:value for header,value in row[1].items()}\
                   for row in networkcrisprdata.iterrows()}
    copy = network.copy()
    nx.set_node_attributes(copy,crisprlabels)
    return copy

def parseRawNetwork(raw_network_filename,totalgenetrees,crisprdata,internal=False,dfonly=False):
    #totalgenetrees = len(os.listdir(newicktreedir))
    df = pd.read_csv(raw_network_filename,
                     delimiter='\t',
                     header=None,
                     names=['raw_score','edge','direction'])
    df = df.loc[df['raw_score'] != 0] #filter edges with 0 score
    df['weight'] = df['raw_score']/totalgenetrees
    df['source'],df['target'],df['has_internal_node'] = zip(*df['edge'].apply(parseEdge))
    df = df.drop(['edge'],axis=1)
    if not internal:
        df = df[df['has_internal_node'] == False]
    if dfonly:
        return df
    network = nx.from_pandas_edgelist(df,'source','target',edge_attr=True)
    return annotateCRISPR(network,crisprdata)

if __name__ == '__main__':
    basedir = sys.argv[1]
    processes = 32
#    pathToLua = os.path.join(basedir,'hide_for_linux/score.lua')
#    pathToHide = os.path.join(basedir,'hide_for_linux/lua')
    genetreefilesdir = os.path.join(basedir,'gene_tree_files/trees/')
    genusname = os.path.basename(os.path.normpath(basedir))
    speciestreefilesdir = os.path.join(basedir,'species_tree_files/species_tree_{}/'.format(genusname))
    treedir, totalgenetrees = fixNewickTrees(speciestreefilesdir,genetreefilesdir,processes)
    print(totalgenetrees)
