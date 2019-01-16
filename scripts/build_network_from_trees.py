import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio import Phylo
import matplotlib.cm as cmx
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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
    if re.search('copy',treestr) != None:
        print('has copy')
    else:
        treestr = re.sub(',',', ',treestr)
        treestr = re.sub('\.1', '',treestr)
        open(outpath,'w').write(treestr)
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

def removeLCA(network):
    nolca = nx.Graph()
    nolcaedges = [(u,v,d) for u,v,d in network.edges(data=True) if (re.search('^LCA\(',u) == None) and (re.search('^LCA\(',v) == None)]
    nolca.add_edges_from(nolcaedges)
    return nolca

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
    print(df.columns)
    network = nx.from_pandas_edgelist(df,'source','target',edge_attr=True)
    for e in network.edges(data=True):
        print(e)
        break
    if not internal:
        network = removeLCA(network)
    return annotateCRISPR(network,crisprdata)

#drawing
def makecolors(cmap):
    scheme = cm = plt.get_cmap(cmap)
    cNorm = colors.Normalize(vmin=0,vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=scheme)
    return scalarMap.to_rgba

def makeEdge(edge,pos,colormap):
    x0,y0 = pos[edge[0]]
    x1,y1 = pos[edge[1]]
    color = colormap(edge[2]['percent_score'])
    newedge = dict(
        type='scatter',
        x = tuple([x0,x1,None]),
        y = tuple([y0,y1,None]),
        hoverinfo='text',
        text = tuple([edge[2]['percent_score']]),
        mode='lines',
        line=dict(
            width=20*edge[2]['percent_score'],
            color='rgb{}'.format(color)
            )
        )
    return newedge

def makeNode(node,pos):
    x,y = pos[node]
    color=lambda x:'red' if re.search('^LCA\(',x) != None else 'blue'
    newnode = dict(type='scatter',
                    x=tuple([x]),
                    y=tuple([y]),
                    text=node,
                    mode='markers',
                    hoverinfo='text',
                    marker=dict(
                        color=color(node),
                        size=30,
                   )
              )
    return newnode

def drawNetwork(network,internal=False,ipy=True,colormap=cmx.viridis):
    if not internal:
        network = removeLCA(network)
    pos = nx.spring_layout(network,iterations=200)
    #draw edges
    edge_trace = []
    for edge in network.edges(data=True):
        edge_trace.append(makeEdge(edge,pos,makecolors(colormap)))
    #draw nodes
    node_trace=[]
    for node in network.nodes(data=True):
        node_trace.append(makeNode(node,pos))
    #figure
    fig = go.Figure(data=edge_trace + node_trace,
                    layout=go.Layout(
                        title='Test HGT Network',
                        titlefont=dict(size=16),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        xaxis=dict(showgrid=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False,
                                   showticklabels=False)
                    )
          )
    if ipy:
        py.iplot(fig, filename='test')
    else:
        py.plot(fig, filename='test')
    return None

def main(speciesdir,genedir,luapath,hidepath):
    #set up dirs
    maindir = 'network_files'
    os.system('mkdir -p {}/all_newick_trees'.format(maindir))
    #convert species tree to newick
    speciestreefile = glob.glob('{}/*run2.t'.format(speciesdir))[0]
    species_tree = getTree(speciestreefile)
    species_tree.root_at_midpoint()
    writeNewickWithOutBrLen(species_tree,'{}/all_newick_trees/species.newick'.format(maindir))
    #convert gene trees to newick
    for treefile in tqdm(glob.iglob('{}/**/*run2.t'.format(genedir)),total=len(os.listdir(genedir))):
        gene_tree = getTree(treefile)
        gene_tree.root_at_midpoint()
        basename = os.path.basename(treefile).split('.')[0]
        nwkfile = '{}/all_newick_trees/{}.newick'.format(maindir,basename)
        writeNewickWithOutBrLen(gene_tree,nwkfile)
    return None

if __name__ == '__main__':
    basedir = '/home/sidreed/thesis_SidReed/hide_test/'
#    pathToLua = os.path.join(basedir,'hide_for_linux/score.lua')
#    pathToHide = os.path.join(basedir,'hide_for_linux/lua')
#    genetreefilesdir = os.path.join(basedir,'gene_tree_files/trees/')
#    speciestreefilesdir = os.path.join(basedir,'species_tree_files/species_tree_ehrlichia/')
#    main(speciestreefilesdir,genetreefilesdir,pathToLua,pathToHide)
    testnet = os.path.join(basedir,'hide/realnetwork.txt')
    totalgenetrees = len(os.listdir('network_files/all_newick_trees'))
    network = parseRawNetwork(testnet,totalgenetrees)
    drawNetwork(network)

