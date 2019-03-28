import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.stats as sst
from functools import reduce
from functools import partial
from itertools import combinations
from collections import defaultdict
from multiprocessing.dummy import Pool as ThreadPool
import networkx as nx

crisprdfpath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'

#Data IO
def isCRISPR(acc,andf):
    internal = ('LCA' in acc)
    if internal:
        return 'internal'
#    anrow = andf[andf['Accession Number'] == acc]
#    iscrispr = anrow['isCRISPR'].values
    iscrispr = andf[andf['Accession Number'] == acc]['isCRISPR'].values
    if iscrispr:
        return 'crispr'
    else:
        return 'non-crispr'

def loaddf(path,andf,filterInternal=True):
    df = pd.read_csv(path,sep='~')
    df['isCRISPR_source'] = df['source'].apply(partial(isCRISPR,andf=andf))
    df['isCRISPR_sink'] = df['sink'].apply(partial(isCRISPR,andf=andf))
    if filterInternal:
        return df[(df['isCRISPR_source'] != 'internal') & (df['isCRISPR_sink'] != 'internal')]
    else:
        return df

def getNodeStatus(df,node):
    if not df[df['source'] == node]['isCRISPR_source'].empty:
        status = df[df['source'] == node]['isCRISPR_source'].iloc[0]
    else:
        status = df[df['sink'] == node]['isCRISPR_sink'].iloc[0]
    return status

def dfToNetwork(df,attrdict=None):
    nodelist = set(df['source']).union(set(df['sink']))
    nodeAttrs = {node:getNodeStatus(df,node) for node in nodelist}
    network = nx.from_pandas_edgelist(df,'source','sink',edge_attr=True)
    nx.set_node_attributes(network,nodeAttrs,'isCRISPR')
    if attrdict != None:
        nx.set_edge_attributes(network,attrdict)
    return network

def loadAll(cdir,cdfpath=crisprdfpath,processes=8):
    andf = pd.DataFrame(json.load(open(cdfpath)))
    pool = ThreadPool(processes)
    ldf = partial(loaddf,andf=andf,filterInternal=True)
    files = [os.path.join(cdir,csv) for csv in os.listdir(cdir)]
    dfs = list(tqdm(pool.imap(ldf,files),total=len(files),desc='reps'))
    nets = [dfToNetwork(df) for df in dfs]
    return nets

def computeEdgeSEs(netlist):
    semlist = {}
    union = lambda x,y: set(x).union(set(y))
    edgeset = list(reduce(union,[nx.edges(net) for net in netlist]))
    allnodes = list(reduce(union,[nx.nodes(net) for net in netlist]))
    for edge in tqdm(list(edgeset),desc='edgeSEM'):
        u,v = edge
        weightlist = []
        for n in netlist:
            try:
                weightlist.append(n[u][v]['weight'])
            except KeyError:
                pass
        semlist.update({tuple([u,v]):{'mean_weight':np.mean(weightlist),'sem_weight':sst.sem(weightlist)}})
    return semlist

def computeModularity(net):
    m = sum(dict(nx.degree(net,weight='weight')).values())
    delta = lambda u,v: 1 if u == v else 0
    q = 0
    for u,v,data in net.edges(data=True):
        ave = data['weight'] - nx.degree(net,u,weight='weight')*nx.degree(net,v,weight='weight')/m
        d = delta(data['isCRISPR_source'],data['isCRISPR_sink'])
        q += ave*d
    return q/m

def diffStats_perBS(netlist,netstat,asdict=False):
    statlist = defaultdict(list)
    for net in tqdm(netlist):
        cnodes = [x for x,y in net.nodes(data=True) if y['isCRISPR']=='crispr']
        ncnodes = [x for x,y in net.nodes(data=True) if y['isCRISPR']=='non-crispr']
        if not asdict:
            cstats = list(dict(netstat(net,cnodes,weight='weight')).values())
            ncstats = list(dict(netstat(net,ncnodes,weight='weight')).values())
        else:
            alls = netstat(net,weight='weight')
            cstats = [alls[x] for x in cnodes]
            ncstats = [alls[x] for x in ncnodes]
        statlist['crispr_mean'].append(np.mean(cstats))
        statlist['crispr_sem'].append(sst.sem(cstats))
        statlist['non-crispr_means'].append(np.mean(ncstats))
        statlist['non-crispr_sems'].append(sst.sem(ncstats))
        mwu = sst.mannwhitneyu(cstats,ncstats)
        statlist['mw_stats'].append(mwu.statistic)
        statlist['mw_pvals'].append(mwu.pvalue)
    return statlist

def diffStats_total(netlist,netstat,asdict=False):
    statdict = {}
    cnodes = list(set([x for net in netlist for x,y in net.nodes(data=True) if y['isCRISPR']=='crispr']))
    ncnodes = list(set([x for net in netlist for x,y in net.nodes(data=True) if y['isCRISPR']=='non-crispr']))
    cstats,ncstats = [],[]
    for net in netlist:
        if not asdict:
            cstats.extend(list(dict(netstat(net,cnodes,weight='weight')).values()))
            ncstats.extend(list(dict(netstat(net,ncnodes,weight='weight')).values()))
        else:
            alls = netstat(net,weight='weight')
            cstats.extend([alls[x] for x in cnodes])
            ncstats.extend([alls[x] for x in ncnodes])
    statdict['crispr_mean'] = np.mean(cstats)
    statdict['crispr_sem'] = sst.sem(cstats)
    statdict['non-crispr_means'] = np.mean(ncstats)
    statdict['non-crispr_sems'] = sst.sem(ncstats)
    mwu = sst.mannwhitneyu(cstats,ncstats)
    statdict['mw_stats'] = mwu.statistic
    statdict['mw_pvals'] = mwu.pvalue
    return statdict

def makeReport(netlist,genusname,perbs=False):
    #get network for plotting
    ses = computeEdgeSEs(netlist)
    maxnodes = max([n.number_of_nodes() for n in netlist])
    basenet = [n for n in netlist if n.number_of_nodes() == maxnodes][0]
    nx.set_edge_attributes(basenet,ses)
    #stat report
    report = {'genus':genusname}
    report['modularity'] = [computeModularity(net) for net in tqdm(netlist)]
    report['assortativity'] = [nx.attribute_assortativity_coefficient(net,'isCRISPR') for net in tqdm(netlist)]
    if perbs:
        report['degree'] = diffStats_perBS(netlist,nx.degree,False)
        report['eg_centrality'] = diffStats_perBS(netlist,nx.eigenvector_centrality_numpy,True)
        report['clustering'] = diffStats_perBS(netlist,nx.clustering,True)
    else:
        report['degree'] = diffStats_total(netlist,nx.degree,False)
        report['eg_centrality'] = diffStats_total(netlist,nx.eigenvector_centrality_numpy,True)
        report['clustering'] = diffStats_total(netlist,nx.clustering,True)
    #report['c_vitaltiy'] = na.diffStats_perBS(netlist,nx.closeness_vitality,True)
    return basenet,report

def writeReport(plotnet,report,genusdir):
    genus = os.path.basename(genusdir)
    pde = nx.to_pandas_edgelist(plotnet)
    pde.to_json(os.path.join(genusdir,'{}_plot_network.json'.format(genus)))
    repath = os.path.join(genusdir,'{}_stat_report.json'.format(genus))
    json.dump(report,open(repath,'w'))
    return None

if __name__ == '__main__':
    genusdir = sys.argv[1]
    perbs = bool(int(sys.argv[2]))
    print(perbs)
    csvs = os.path.join(genusdir,'network_files/csvs')
    print(csvs)
    genus = os.path.basename(genusdir)
    print(genus)
    netlist = loadAll(csvs)
    plotnet, report = makeReport(netlist,genusdir,perbs)
    writeReport(plotnet,report,genusdir)

#Node stats (compares crispr vs non-crispr with test)
#edge weight
#katz centrality
#eigenvector entrality nx.eigenvector_centrality_numpy(net,weight='weight')
#closeness vitality nx.closeness_vitality(net,weight='weight')
#clustering
#Ntwork stats
#modularity computeModularity(df)
#assortativity attribute_assortativity_coefficient(net,'isCRISPR')
#dfcols are
#source, sink, raw_weight, weight, direction, isCRISPR_source, isCRISPR_sink

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
#def computeModularity(df):
#    m = sum(df['weight']) #total edge weight
#    delta = lambda u,v: 1 if u == v else 0
#    tipnodes = getTipNodeList(df)
#    q = 0
#    for (u,v) in combinations(tipnodes,2):
#        ave = getEdgeWeight(df,u,v) - getNodeDeg(df,u)*getNodeDeg(df,v)/(2*m)
#        d = delta(getNodeStatus(df,u),getNodeStatus(df,u))
#        q += ave*d
#    return q/(2*m)
#
#
#def getTipNodeList(df):
#    sourcetips = df[df['isCRISPR_source'] != 'internal']['source']
#    sinktips = df[df['isCRISPR_sink'] != 'internal']['sink']
#    tipnodes = list(set(sourcetips).union(set(sinktips)))
#    return tipnodes
#
#def edgesByStatus(df,status):
#    return df[(df['isCRISPR_source'] == status) | (df['isCRISPR_sink'] == status)]
#
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
#
##Node Functions
#def getNodeDf(df,nodename):
#    return df[(df['source'] == nodename) | (df['sink'] == nodename)]
#
#def getNodeDeg(df,nodename):
#    return sum(getNodeDf(df,nodename)['weight'])
#
##Edge Functions
#def getEdgeDf(dfconcat,u,v):
#    forward = df[(df['source'] == u) & (df['sink'] == v)]
#    reverse = df[(df['source'] == v) & (df['sink'] == u)]
#    return pd.concat([forward,reverse])
#
#def getEdgeWeight(df,u,v):
#    edge = df[(df['source'] == u) & (df['sink'] == v)]
#    if not edge.empty:
#        return edge['weight'].values[0]
#    else:
#        return 0
#
##Statistics
#def avgEdgeWeight(df):
#    return sum(df['weight'])/df.size[0]
#
#def computeModularity(net):
#    m = sum(dict(nx.number_of_edges(net,weight='weight')))
#    delta = lambda u,v: 1 if u == v else 0
#    q = 0
#    for u,v,data in net.edges(data=True):
#        ave = data['weight'] - nx.degree(net,u,weight='weight')*nx.degree(net,u,weight='weight')
#        d = delta(data['isCRISPR_source'],data['isCRISPR_sink'])
#        q += ave*d
#    return q/(m)
#
#def loadAll(cdir,cdfpath=crisprdfpath,processes=32,asnet=False):
#    andf = pd.DataFrame(json.load(open(cdfpath)))
#    pool = ThreadPool(processes)
#    ldf = partial(loaddf,andf=andf)
#    files = [os.path.join(cdir,csv) for csv in os.listdir(cdir)]
#    dfs = list(tqdm(pool.imap(ldf,files),total=len(files),desc='dfs'))
#    return dfs

