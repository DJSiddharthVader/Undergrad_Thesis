import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import networkx as nx
from Bio import Phylo
from tqdm import tqdm

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
def getTree(treefilename):
    return list(Phylo.parse(treefilename,'nexus'))[-1]

def writeNewickWithOutBrLen(treeobj,outpath):
    treestr = re.sub(':\d\.\d+','',treeobj.format('newick'))
    hascopy = 0
    if re.search('copy',treestr) != None:
        hascopy = 1
    else:
        treestr = re.sub(',',', ',treestr)
        treestr = re.sub('\.1', '',treestr)
        open(outpath,'w').write(treestr)
    return hascopy

def fixNewickTrees(speciesdir,genedir):
    #set up dirs
    maindir = 'network_files'
    os.system('mkdir -p {}/all_newick_trees'.format(maindir))
    #convert species tree to newick
    speciestreefile = glob.glob('{}/*run2.t'.format(speciesdir))[0]
    species_tree = getTree(speciestreefile)
    species_tree.root_at_midpoint()
    writeNewickWithOutBrLen(species_tree,'{}/all_newick_trees/species.newick'.format(maindir))
    #convert gene trees to newick
    hascopy = 0
    for treefile in tqdm(glob.iglob('{}/**/*run2.t'.format(genedir)),total=len(os.listdir(genedir))):
        gene_tree = getTree(treefile)
        gene_tree.root_at_midpoint()
        basename = os.path.basename(treefile).split('.')[0]
        nwkfile = '{}/all_newick_trees/{}.newick'.format(maindir,basename)
        hascopy += writeNewickWithOutBrLen(gene_tree,nwkfile)
    print('{} of {} gene trees had duplicate taxa, not considered'.format(hascopy,len(os.listdir(genedir))))
    numgenetrees= len(os.listdir('network_files/all_newick_trees')) -1 #-1 for the species tree
    return '{}/all_newick_trees'.format(maindir), numgenetrees

def runHideOnTreeDir(treedir,hiderDir):
    luaexe = os.path.join(hideDir,'lua')
    networktxt = subprocess.run(['mb'],
                            stdin=mbscript,
                            stdout=logfile,
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

def parseRawNetwork(raw_network_filename,totalgenetrees,crisprdata,internal=False):
    #totalgenetrees = len(os.listdir(newicktreedir))
    df = pd.read_csv(raw_network_filename,
                     delimiter='\t',
                     header=None,
                     names=['raw_score','edge','direction'])
    df = df.loc[df['raw_score'] != 0] #filter edges with 0 score
    df['percent_score'] = df['raw_score']/totalgenetrees
    df['source'],df['target'],df['has_internal_node'] = zip(*df['edge'].apply(parseEdge))
    df = df.drop(['edge'],axis=1)
    if not internal:
        df = df[df['has_internal_node'] == False]
    network = nx.from_pandas_edgelist(df,'source','target',edge_attr=True)
    return annotateCRISPR(network,crisprdata)

if __name__ == '__main__':
#    basedir = '/home/sidreed/thesis_SidReed/hide_test/'
    basedir = sys.argv[1]
#    pathToLua = os.path.join(basedir,'hide_for_linux/score.lua')
#    pathToHide = os.path.join(basedir,'hide_for_linux/lua')
    genetreefilesdir = os.path.join(basedir,'gene_tree_files/trees/')
    genusname = os.path.basename(os.path.normpath(basedir))
    speciestreefilesdir = os.path.join(basedir,'species_tree_files/species_tree_{}/'.format(genusname))
    treedir, totalgenetrees = fixNewickTrees(speciestreefilesdir,genetreefilesdir)
    print(totalgenetrees)
#    testnet = runHideOnTreeDir(treedir)
#    network = parseRawNetwork(testnet,totalgenetrees)
