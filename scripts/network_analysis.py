import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
import scipy.stats as sst
from itertools import combinations
from functools import reduce,partial
from collections import defaultdict,Counter
from multiprocessing.dummy import Pool as ThreadPool

#crisprdfpath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'
crisprdfpath = '/home/sid/thesis_SidReed/data/all_data.tsv'
toolonggenera = ['Bifidobacterium',
                 'Chlamydia',
                 'Legionella',
                 'Neisseria',
                 'Streptomyces',
                 'Yersinia'] #genera that take too long to calculate closeness_vitality, skip it for those genera
#Node stats (compares crispr vs non-crispr with test)
#edge weight
#katz centrality
#eigenvector entrality nx.eigenvector_centrality_numpy(net,weight='weight')
#closeness vitality nx.closeness_vitality(net,weight='weight')
#clustering
#Ntwork stats
#modularity modularity(df)
#assortativity attribute_assortativity_coefficient(net,'isCRISPR')
#dfcols are
#source, sink, raw_weight, weight, direction, isCRISPR_source, isCRISPR_sink


#Utilities
def isCRISPR(acc,andf):
    #determine whether an OTU (accession number) has a crispr system
    internal = ('LCA' in acc)
    if internal:
        return 'internal'
    iscrispr = andf[andf['Accession Number'] == acc]['isCRISPR'].values
    if (iscrispr.size > 0):
        if iscrispr.all():
            return 'crispr'
        else:
            return 'non-crispr'
    else:
        return 'non-crispr'

def getNodeStatus(df,node):
    #determine if the node in an edgelist is crispr or not
    if not df[df['source'] == node]['isCRISPR_source'].empty:
        status = df[df['source'] == node]['isCRISPR_source'].iloc[0]
    else:
        status = df[df['sink'] == node]['isCRISPR_sink'].iloc[0]
    return status

#Data IO
def loaddf(path,andf,filterInternal=True):
    #load the edgelists produced by HiDe and filter out the edges connecting at least 1 internal (non-leaf) node
    df = pd.read_csv(path,sep='~')
    if df.shape[0] == 0:
        return False
    df['isCRISPR_source'] = df['source'].apply(partial(isCRISPR,andf=andf))
    df['isCRISPR_sink'] = df['sink'].apply(partial(isCRISPR,andf=andf))
    if filterInternal:
        return df[(df['isCRISPR_source'] != 'internal') & (df['isCRISPR_sink'] != 'internal')]
    else:
        return df

def dfToNetwork(df,attrdict=None):
    #convert a dataframe of edgelists from HiDe into a networkx network object
    nodelist = set(df['source']).union(set(df['sink']))
    nodeAttrs = {node:getNodeStatus(df,node) for node in nodelist}
    network = nx.from_pandas_edgelist(df,'source','sink',edge_attr=True)
    nx.set_node_attributes(network,nodeAttrs,'isCRISPR')
    if attrdict != None:
        nx.set_edge_attributes(network,attrdict)
    return network

def loadNetwork(networkfile,andf,filterInternal):
    netdf = loaddf(networkfile,andf,filterInternal)
    network = dfToNetwork(netdf)
    return network

def loadAllNetworks(cdir,processes,genus):
    #load all the HiDe edgelists and convert them to networkx objects and return the list of networks
    if not(os.path.exists(cdir)):
        print('no network csv dir')
        sys.exit()
    files = [os.path.join(cdir,csv) for csv in os.listdir(cdir)]
    andf = pd.read_csv(crisprdfpath,sep='\t')
    ln = partial(loadNetwork,andf=andf,filterInternal=True)
    if (len(files)) == 0:
        print('No network files')
        sys.exit()
    pool = ThreadPool(processes)
    nets = list(tqdm(pool.imap(ln,files),total=len(files),desc='loading networks...'))
    nets = [net for net in nets if isinstance(net,pd.core.frame.DataFrame)]
    return nets

#Stats
def computeEdgeSEs(netlist):
    #for a list of networks compute the average weight and sd. err. for each edge (pair of nodes)
    semlist = {}
    union = lambda x,y: set(x).union(set(y))
    edgeset = list(reduce(union,[nx.edges(net) for net in netlist]))
    allnodes = list(reduce(union,[nx.nodes(net) for net in netlist]))
    for edge in tqdm(list(edgeset),desc='edge Mean/SE'):
        u,v = edge
        weightlist = []
        for n in netlist:
            try:
                weightlist.append(n[u][v]['weight'])
            except KeyError:
                pass
        semlist.update({tuple([u,v]):{'mean_weight':np.mean(weightlist),'sem_weight':sst.sem(weightlist)}})
        #for each edge calculate the average and std. err. of the edge weight across the 1000 replicates
    return semlist

def modularity(net):
    #the formula is here https://en.wikipedia.org/wiki/Louvain_modularity
    #modularity Q = 1/2m \sum{all_edges} [edge - degree(u)*degree(v)/2(m)]*1 \delta(u,v)
    #where u,v are the nodes of an edge, m is the sum of all edges in the network and
    #\delta(u,v) is the knonecker delta if u and v are both crispr and non-crispr
    #Q is between -1 and 1, with 1 indicating high community structure (crispr nodes are more strongly connected to crispr nodes, same for non-crispr nodes)
    q = 0
    m = sum(dict(nx.degree(net,weight='weight')).values())
    delta = lambda u,v: 1 if u == v else 0
    for u,v,data in net.edges(data=True): #iterate over all edges
        ave = data['weight'] - nx.degree(net,u,weight='weight')*nx.degree(net,v,weight='weight')/(2*m)
        d = delta(data['isCRISPR_source'],data['isCRISPR_sink'])
        q += ave*d
    return q/(2*m)

def computeDensity(net):
    #each edge has a max value of 1 since values were normalized during network construction
    #density is the sum of all edges or the maximum possible sum of all edges (1*number of possible edges)
    m = sum(dict(nx.degree(net,weight='weight')).values()) #sum of all edges in the graph
    n = nx.number_of_nodes(net)*(nx.number_of_nodes(net)-1) #number of possible edges (i.e max possible total edge weight)
    return m/n

def computeStat(network,statfnc,fmt):
    cnodes = list(set([x for x,y in network.nodes(data=True) if y['isCRISPR']=='crispr']))
    ncnodes = list(set([x for x,y in network.nodes(data=True) if y['isCRISPR']=='non-crispr']))
    if fmt == 'edict':
        crispr = list(dict(statfnc(network,cnodes,weight='weight')).values())
        noncrispr = list(dict(statfnc(network,ncnodes,weight='weight')).values())
    elif fmt == 'ndict':
        stats =  dict(statfnc(network,weight='weight'))
        crispr = list({k:v for k,v in stats.items() if k in cnodes}.values())
        noncrispr = list({k:v for k,v in stats.items() if k in ncnodes}.values())
    elif fmt == 'list':
        alls = statfnc(network,weight='weight') #a statistic for each node in the network
        crispr = [alls[x] for x in cnodes]
        noncrispr = [alls[x] for x in ncnodes]
    elif fmt == 'subgraph':
        cnet = nx.subgraph(network,nbunch=cnodes)   #take the subgraph of only crispr nodes
        ncnet = nx.subgraph(network,nbunch=ncnodes) #take the subgraph of only non-crispr nodes
        try:
            crispr = cnet.size(weight='weight')/cnet.number_of_edges()
        except ZeroDivisionError:
            crispr = 0 #only reached if no crispr nodes in the network
        try:
            noncrispr = ncnet.size(weight='weight')/ncnet.number_of_edges()
        except ZeroDivisionError:
            noncrispr = 0 #only reached if no non-crispr nodes in the network
    else:
        raise ValueError('invalid calculation')
    return {'c':crispr,'nc':noncrispr}

def separateStats(statpairs):
    if isinstance(list(statpairs[0].values())[0],list):
        cstats = [y for x in statpairs for y in x['c']]
        ncstats = [y for x in statpairs for y in x['nc']]
        return (cstats,ncstats)
    elif isinstance(list(statpairs[0].values())[0],(int,float,bool)):
        cstats = [x['c'] for x in statpairs]
        ncstats = [x['nc'] for x in statpairs]
        return (cstats,ncstats)
    else:
        raise ValueError('invalid statpairs')

def diffStats_total(netlist,statfnc,fmt,processes):
    #here statfnc is a function (mostly from the networkx package) that will return a statistic (single value, list, dict) given a single network
    #Then seperate the calculated statistics into only the crispr and non-crispr nodes
    #finally reutrn the list of statistics
    #   - the crispr mean and SE
    #   - the noncrispr mean,SE and
    #   - the Mann-Whitney U stat,pval testing if the the crispr and non-crisprs lists are different
    pool = ThreadPool(processes)
    calcstat = partial(computeStat,statfnc=statfnc,fmt=fmt)
    statpairs = list(tqdm(pool.imap(calcstat,netlist),total=len(netlist),desc=statfnc.__name__))
    cstats, ncstats = separateStats(statpairs)
    try:
        mwu = sst.mannwhitneyu(cstats,ncstats)
        return {'crispr_mean':np.mean(cstats),
                'crispr_sem':sst.sem(cstats),
                'non-crispr_mean':np.mean(ncstats),
                'non-crispr_sem':sst.sem(ncstats),
                'mw_stat':mwu.statistic,
                'mw_pval':mwu.pvalue
               }
    except ValueError: #All numbers are identical in mannwhitneyu
        return {'crispr_mean':np.mean(cstats),
                'crispr_sem':sst.sem(cstats),
                'non-crispr_mean':np.mean(ncstats),
                'non-crispr_sem':sst.sem(ncstats),
                'mw_stat':-1,
                'mw_pval':-1
               }

def makeReport(netlist,genus,processes):
    ses = computeEdgeSEs(netlist) #get the mean and SE for all edge weights
    maxnodes = max([n.number_of_nodes() for n in netlist])
    basenet = [n for n in netlist if n.number_of_nodes() == maxnodes][0] #network with the most nodes, i.e. the one with all OTUs in the genus in it
    nx.set_edge_attributes(basenet,ses) #set the network to have the edge values be the mean across all bootstraps
    #Create a dictionary of statistics for the set of networks
    report = {'genus':genus}
    report['degree'] = diffStats_total(netlist,nx.degree,'edict',processes)
    report['clustering'] = diffStats_total(netlist,nx.clustering,'list',processes)
    report['mean_edge_weight'] = diffStats_total(netlist,nx.degree,'subgraph',processes)
    if genus in toolonggenera:
        report['closeness_vitaltiy'] = {'crispr_mean':'takes_too_long', 'crispr_sem':'takes_too_long', 'non-crispr_mean':'takes_too_long', 'non-crispr_sem':'takes_too_long', 'mw_stat':'takes_too_long', 'mw_pval':'takes_too_long'}
    else:
        report['closeness_vitaltiy'] = diffStats_total(netlist,nx.closeness_vitality,'ndict',processes)
    report['centrality'] = diffStats_total(netlist,nx.eigenvector_centrality_numpy,'ndict',processes)
    report['assortativity'] = [nx.attribute_assortativity_coefficient(net,'isCRISPR') \
                               for net in tqdm(netlist,desc='assortativity')]
    report['modularity'] = [modularity(net) for net in tqdm(netlist,desc='modularity')]
    return basenet,report

def main(genusdir,sttype,processes):
    genus = os.path.basename(genusdir)
    plotfile = './{}_{}_plot_network.json'.format(genus,sttype)
    statfile = './{}_{}_stat_report.json'.format(genus,sttype)
    if os.path.isfile(plotfile) and os.path.isfile(statfile):
        print('already calculated')
        return None
    csvs = os.path.join(genusdir,'network_files_{}/csvs'.format(sttype))
    netlist = loadAllNetworks(csvs,processes,genus)
    plotnet, report = makeReport(netlist,genus,processes)
    pde = nx.to_pandas_edgelist(plotnet)
    pde.to_json(plotfile) #writes the average network to json object so it can be loaded/plotted
    json.dump(report,open(statfile,'w')) #writes all of the computed statistics to a json file
    return None


if __name__ == '__main__':
    if len(sys.argv) > 1:
        sttype = sys.argv[1] #WGS or 16S
    else:
        sttype = 'WGS'
    if len(sys.argv) > 2:
        genusdir = os.path.abspath(sys.argv[2])
    else:
        genusdir = os.getcwd()
    if len(sys.argv) > 3:
        processes = int(sys.argv[3])
    else:
        processes = 48
    print(genusdir,sttype,processes)
    if sttype not in ['16S','WGS']:
        print('invalid arg 1, use {16S|WGS}')
        sys.exit()
    main(genusdir,sttype,processes)


#DEPRECIARED
#def avgEdgeWeight(netlist):
#    alledges = pd.concat(netlist)
#    cweights = alledges[(alledges['isCRISPR_source'] == 'crispr') | (alledges['isCRISPR_sink'] == 'crispr')]['weight']
#    ncweights = alledges[~((alledges['isCRISPR_source'] == 'crispr') | (alledges['isCRISPR_sink'] == 'crispr' ))]['weight']
#    cmean = np.mean(cweights)
#    cse = sst.sem(cweights)
#    ncmean = np.mean(ncweights)
#    ncse = sst.sem(ncweights)
#    ustat = sst.mannwhitneyu(cweights,ncweights,alternative='two-sided')
#    return {'crispr_mean_edge_weight':cmean,
#            'non_crispr_mean_edge_weight':ncmean,
#            'crispr_std_err':cse,
#            'non_crispr_std_err':ncse,
#            'mann_whitney_u_stat':ustat.statistic,
#            'mann_whitney_u_pvalue':ustat.pvalue}
#def modularity(df):
#    m = sum(df['weight']) #total edge weight
#    delta = lambda u,v: 1 if u == v else 0
#    tipnodes = getTipNodeList(df)
#    q = 0
#    for (u,v) in combinations(tipnodes,2):
#        ave = getEdgeWeight(df,u,v) - getNodeDeg(df,u)*getNodeDeg(df,v)/(2*m)
#        d = delta(getNodeStatus(df,u),getNodeStatus(df,u))
#        q += ave*d
#    return q/(2*m)
#def getTipNodeList(df):
#    sourcetips = df[df['isCRISPR_source'] != 'internal']['source']
#    sinktips = df[df['isCRISPR_sink'] != 'internal']['sink']
#    tipnodes = list(set(sourcetips).union(set(sinktips)))
#    return tipnodes
#def edgesByStatus(df,status):
#    return df[(df['isCRISPR_source'] == status) | (df['isCRISPR_sink'] == status)]
#def computeEdgeSEs(netlist):
#    semlist = {}
#    union = lambda x,y: set(x).union(set(y))
#    edgeset = list(reduce(union,[nx.edges(net) for net in netlist]))
#    allnodes = list(reduce(union,[nx.nodes(net) for net in netlist]))
#    for edge in list(edgeset):
#        u,v = edge
#        weightlist = []
#        for n in netlist:
#            try:
#                weightlist.append(n[u][v]['weight'])
#            except KeyError:
#                pass
#        semlist.update({tuple([u,v]):{'mean_weight':np.mean(weightlist),'sem_weight':sst.sem(weightlist)}})
#    return semlist
#def getNodeDf(df,nodename):
#    return df[(df['source'] == nodename) | (df['sink'] == nodename)]
#def getNodeDeg(df,nodename):
#    return sum(getNodeDf(df,nodename)['weight'])
#def getEdgeDf(dfconcat,u,v):
#    forward = df[(df['source'] == u) & (df['sink'] == v)]
#    reverse = df[(df['source'] == v) & (df['sink'] == u)]
#    return pd.concat([forward,reverse])
#def getEdgeWeight(df,u,v):
#    edge = df[(df['source'] == u) & (df['sink'] == v)]
#    if not edge.empty:
#        return edge['weight'].values[0]
#    else:
#        return 0
#def avgEdgeWeight(df):
#    return sum(df['weight'])/df.size[0]
#def modularity(net):
#    m = sum(dict(nx.number_of_edges(net,weight='weight')))
#    delta = lambda u,v: 1 if u == v else 0
#    q = 0
#    for u,v,data in net.edges(data=True):
#        ave = data['weight'] - nx.degree(net,u,weight='weight')*nx.degree(net,u,weight='weight')
#        d = delta(data['isCRISPR_source'],data['isCRISPR_sink'])
#        q += ave*d
#    return q/(m)
#def loadAllNetworks(cdir,processes,cdfpath=crisprdfpath):
#    #load all the HiDe edgelists and convert them to network x objects and return the list of networks
#    #andf = pd.DataFrame(json.load(open(cdfpath)))
#    andf = pd.read_csv(cdfpath,sep='\t')
#    pool = ThreadPool(processes)
#    ldf = partial(loaddf,andf=andf,filterInternal=True)
#    files = [os.path.join(cdir,csv) for csv in os.listdir(cdir)]
#    dfs = list(tqdm(pool.imap(ldf,files),total=len(files),desc='loading HiDes...'))
#    nets = [dfToNetwork(df) for df in dfs]
#    return nets
#def modularity(net):
#    q = 0
#    m = sum(dict(nx.degree(net,weight='weight')).values())
#    delta = lambda u,v: 1 if u == v else 0
#    for u,v,data in net.edges(data=True):
#        ave = data['weight'] - nx.degree(net,u,weight='weight')*nx.degree(net,v,weight='weight')/m
#        d = delta(data['isCRISPR_source'],data['isCRISPR_sink'])
#        q += ave*d
#    return q/m
#def diffStats_perBS(netlist,netstat,fmt=False):
#    statlist = defaultdict(list)
#    for net in tqdm(netlist):
#        cnodes = [x for x,y in net.nodes(data=True) if y['isCRISPR']=='crispr']
#        ncnodes = [x for x,y in net.nodes(data=True) if y['isCRISPR']=='non-crispr']
#        if fmt == 'dict':
#            cstats = list(dict(netstat(net,cnodes,weight='weight')).values())
#            ncstats = list(dict(netstat(net,ncnodes,weight='weight')).values())
#        elif fmt == 'list':
#            alls = netstat(net,weight='weight')
#            cstats = [alls[x] for x in cnodes]
#            ncstats = [alls[x] for x in ncnodes]
#        elif fmt == 'subgraph':
#            cnet = nx.edges(net,nbunch=cnodes)
#            cstats.append(cnet.size(weight='weight')/cnet.number_of_edges())
#            ncnet = nx.edges(net,nbunch=ncnodes)
#            cstats.append(ncnet.size(weight='weight')/ncnet.number_of_edges())
#        else:
#            raise ValueError('invalid calculation')
#    mwu = sst.mannwhitneyu(cstats,ncstats)
#    return {'crispr_mean':np.mean(cstats)
#            'crispr_sem':sst.sem(cstats)
#            'non-crispr_mean':np.mean(ncstats)
#            'non-crispr_sem':sst.sem(ncstats)
#            'mw_stat':mwu.statistic
#            'mw_pval':mwu.pvalue
#           }
#def diffStats_total(netlist,netstat,fmt='edict'):
#    #here netstat is a function (mostly from the networkx package) that will return a statistic (single value, list, dict) given a single network
#    #Then seperate the calculated statistics into only the crispr and non-crispr nodes
#    #finally reutrn the list of statistics
#    #   - the crispr mean and SE
#    #   - the noncrispr mean,SE and
#    #   - the Mann-Whitney U stat,pval testing if the the crispr and non-crisprs lists are different
#    cstats,ncstats = [],[]
#    for net in tqdm(netlist):
#        cnodes = list(set([x for x,y in net.nodes(data=True) if y['isCRISPR']=='crispr']))
#        ncnodes = list(set([x for x,y in net.nodes(data=True) if y['isCRISPR']=='non-crispr']))
#        if fmt == 'edict':
#            cstats.extend(list(dict(netstat(net,cnodes,weight='weight')).values()))
#            ncstats.extend(list(dict(netstat(net,ncnodes,weight='weight')).values()))
#        elif fmt == 'ndict':
#            stats =  dict(netstat(net,weight='weight'))
#            cstats.extend(list({k:v for k,v in stats.items() if k in cnodes}.values()))
#            ncstats.extend(list({k:v for k,v in stats.items() if k in ncnodes}.values()))
#        elif fmt == 'list':
#            alls = netstat(net,weight='weight') #a statistic for each node in the network
#            cstats.extend([alls[x] for x in cnodes])
#            ncstats.extend([alls[x] for x in ncnodes])
#        elif fmt == 'subgraph':
#            cnet = nx.subgraph(net,nbunch=cnodes)   #take the subgraph of only crispr nodes
#            ncnet = nx.subgraph(net,nbunch=ncnodes) #take the subgraph of only non-crispr nodes
#            try:
#                cstats.append(cnet.size(weight='weight')/cnet.number_of_edges())
#            except ZeroDivisionError:
#                cstats.append(0) #only reached if no crispr nodes in the network
#            try:
#                ncstats.append(ncnet.size(weight='weight')/ncnet.number_of_edges())
#            except ZeroDivisionError:
#                ncstats.append(0) #only reached if no non-crispr nodes in the network
#        else:
#            raise ValueError('invalid calculation')
#    mwu = sst.mannwhitneyu(cstats,ncstats)
#    return {'crispr_mean':np.mean(cstats),
#            'crispr_sem':sst.sem(cstats),
#            'non-crispr_mean':np.mean(ncstats),
#            'non-crispr_sem':sst.sem(ncstats),
#            'mw_stat':mwu.statistic,
#            'mw_pval':mwu.pvalue
#           }

