import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.stats as sst
from functools import reduce

missing_val = -1
marko_filename = 'markophylo_results.txt'
marko_file_errors = ['[1] "plotting tree"','"Error']
report_filename = 'markophylo_results.txt'

#Output Description
"""
For each statistic there are 6 values returned
- the mean of the statistic of all CRISPR taxa in the genus across all bootstraps
- the std. error of the statistic of all CRISPR taxa in the genus across all bootstraps
- the mean of the statistic of all Non-CRISPR taxa in the genus across all bootstraps
- the std. error of the statistic of all Non-CRISPR taxa in the genus across all bootstraps
- the U statistic from the Mann-Whiteney U test between the CRISPR and Non-CRISPR values
-  the p-value from the Mann-Whiteney U test between the CRISPR and Non-CRISPR values
*Note: For statistics on the entire network (modularity, assortativity) only the mean
       and std. error are computed as there is no comparison to be made
A list of the statistics and descriptions (nx = networkx function, partitions are CRISPR/Non-CRISPR)
## Indel Rate
- gene insertion/deletion rates for a partition of nodes in a species tree
- partitons are CRISPR, Non-CRISPR and internal (non-leaf nodes in the tree, unused, avoids confounding results)
- calculated using markophylo R package using a 16S tree as the species tree and a gene P/A matrix
## Degree
- degree is the sum of all connected edge weights of a single node
- mean of degree of all nodes in a partition across all bootstraps
- degree includes connections between CRISPR and Non-CRISPR nodes
- calculated using the nx.degree function applied to the entire network
## Mean Edge Weight
 - compute the mean edge weight for a partition of nodes across all boostraps
 - calculated using the nx.degree fuction applied to CRISPR/Non-CRISPR subgraphs
## Weighted Clustering Coefficient
 - clustering coefficient of a node is weighted fraction of all triangles that node could exist in
 - take the mean of the WCC of all nodes in a partition across all bootstraps
 - intended to measure how tightly connected nodes (how much/how exclusive) are to each other
 - calculated using the nx.clustering fuction applied to each node
 - [link](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.cluster.clustering.html#networkx.algorithms.cluster.clustering)
## Closeness Vitality
 - the sum of the change in distance between all nodes if a given node is excluded
 - how much does a taxa contribute to transfer across the network (kinda loose)
 - calculated using nx.closeness_vitality for each node
 - [link](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.vitality.closeness_vitality.html)
## Eigenvector Centrality
 - a weighted version of degree centraility where the weight of each node is it's degree centrality
 - are there "hub" taxa? or sets of highly connected taxa? or is it relatively  even
 - calculated using nx.eigenvector_centrality_numpy for each node
 - [link](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.eigenvector_centrality_numpy.html#networkx.algorithms.centrality.eigenvector_centrality_numpy)
## Modularity
 - given 2 edge partitions, the relative density of edges within a partition to those outside of it
 - are CRISPR and Non-CRISPR taxa distinct "communities" based on their edges? or does trasnfer occur between all taxa (mixing)?
 - calculated using my implementation of modularity, based on [this definition](https://en.wikipedia.org/wiki/Louvain_modularity)
## Assortativity
 - "Assortativity measures the similarity of connections in the graph with respect to the given attribute."
 - are the connections made by 2 CRISPR/Non-CRISPR nodes different?
 - calculated using nx.attribute_assortativity_coefficient for each node
 - [link](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.assortativity.attribute_assortativity_coefficient.html#networkx.algorithms.assortativity.attribute_assortativity_coefficient)
"""


def hasRate(lines,cat):
    boollist = [True if cat in l else False for l in lines ]
    return reduce(lambda x,y: x or y, boollist)

def getRates(lines,pattern,has_c,has_nc):
    if pattern not in lines:
        return missing_val,missing_val,missing_val
    rates = lines[lines.index(pattern)+4].split(' ')[1:]
    if set(marko_file_errors).intersection(set(rates)):
        return missing_val,missing_val,missing_val
    rates = [float(x) for x in rates if x != '']
    if has_c and has_nc:
        crate,ncrate,intrate = rates
    elif has_c:
        crate,intrate = rates
        ncrate = missing_val
    elif has_nc:
        ncrate,intrate = rates
        crate = missing_val
    else:
        raise ValueError('Only has internal partition?')
    return crate,ncrate,intrate

def parseMarko(basedir,genus):
    filepath = os.path.join(basedir,genus,marko_filename)
    if not(os.path.exists(filepath)):
        return None
    lines = [x.strip() for x in open(filepath).readlines()]
    if len(lines) == 0:
        return None
    has_c = hasRate(lines,'$crispr')
    has_nc = hasRate(lines,'$non_crispr')
    crate,ncrate,intrate = getRates(lines,'$rates',has_c,has_nc)
    csem,ncsem,intsem = getRates(lines,'$se$rates',has_c,has_nc)
    return {'genus':genus,
            'indel_rate.crispr_rate':crate,
            'indel_rate.non.crispr_rate':ncrate,
            'indel_rate.internal_rate':intrate,
            'indel_rate.crispr_sem':csem,
            'indel_rate.non.crispr_sem':ncsem,
            'indel_rate.internal_sem':intsem}

def makeMarkoDF(basedir,inputs):
    marko_data=[parseMarko(basedir,genus) for genus in tqdm(inputs)]
    marko_df = pd.DataFrame(list(filter(None,marko_data)))
    return marko_df

def reportToDF(report_file):
    if not(os.path.exists(report_file)):
        return pd.DataFrame()
    report = json.load(open(report_file))
    report['modularity.mean'] = np.mean(report['modularity'])
    report['modularity.sem'] = sst.sem(report['modularity'])
    report['assortativity.mean'] = np.mean(report['assortativity'])
    report['assortativity.sem'] = sst.sem(report['assortativity'])
    report.pop('modularity')
    report.pop('assortativity')
    return pd.io.json.json_normalize(report)

def makeReportDF(basedir,genera,smode):
    filepaths = [os.path.join(basedir,genus,'{}_{}_stat_report.json'.format(genus,smode)) for genus in genera]
    dfs = [reportToDF(report_file) for report_file in tqdm(filepaths)]
    report_df = pd.concat([x for x in dfs if not(x.empty)])
    return report_df

def main(basedir,marko_outpath,report_outpath,smode):
    inputs = [genus for genus in os.listdir(basedir) \
              if os.path.isdir(os.path.join(basedir,genus))]
    marko_df = makeMarkoDF(basedir,inputs)
    marko_df.to_csv(marko_outpath,index=False)
    report_df = makeReportDF(basedir,inputs,smode)
    report_df.to_csv(report_outpath,index=False)
    combined_df = marko_df.merge(report_df,how='outer',on='genus')
    combined_df.to_csv('all_{}_results.csv'.format(smode),index=False)
    return None


if __name__ == '__main__':
    if len(sys.argv) > 1:
        basedir = sys.argv[1]
    else:
        basedir =  os.getcwd()
    if len(sys.argv) > 2:
        marko_outpath = sys.argv[2]
    else:
        marko_outpath = 'all_markophylo_results.csv'
    if len(sys.argv) > 3:
        smode = sys.argv[3]
    else:
        smode = 'WGS'
    if len(sys.argv) > 4:
        report_outpath= sys.argv[4]
    else:
        report_outpath = 'all_{}_report_results.csv'.format(smode)
    print(basedir,marko_outpath,smode,report_outpath)
    main(basedir,marko_outpath,report_outpath,smode)

