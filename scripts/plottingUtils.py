import os
import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from functools import partial

def loadnetwork(path):
    netjson = json.load(open(path))
    df = pd.DataFrame(netjson)
    df = df[df['source'] != df['sink']].drop('target',axis=1)
    net = nx.from_pandas_edgelist(df,'source','sink',edge_attr=True)
    return df,net

def getNodeStatus(df):
    sos = {s:list(set(df[df['source'] == s]['isCRISPR_source']))[0] for s in set(df['source'])}
    sodf = pd.DataFrame.from_dict(sos,orient='index')
    sodf.columns = ['source']
    sis = {s:list(set(df[df['sink'] == s]['isCRISPR_sink']))[0] for s in set(df['sink'])}
    sidf = pd.DataFrame.from_dict(sis,orient='index')
    sidf.columns = ['sink']
    cdf = pd.concat([sodf,sidf],axis=1,sort=False).fillna(0)
    consensus = lambda row: row['source'] if (row['sink'] == 0) else row['sink']
    cdf['status'] =  cdf.apply(consensus,axis=1)
    return {i:row['status'] for i,row in cdf.iterrows()}

def plotNetwork(net,df,width_scaling=4000,alpha=0.8,legend_entries=6,legend_decimals=6,ccol='blue',nccol='red',edge_cmap=plt.cm.viridis_r,layout=nx.spring_layout,seed=9,entry_width=0.03*0.63,left_legend_offset=1.05,save=False):
    #set up plot
    np.random.seed(seed)
    fig,ax = plt.subplots(1,1,figsize=(20,15))
    ax.set_axis_off()
    nodeidx = {node:i for i,node in enumerate(net.nodes())}
    #color nodes
    nodestatus = getNodeStatus(df)
    cnodes = [n for n in net.nodes() if nodestatus[n] == 'crispr' ]
    ncnodes = [n for n in net.nodes() if nodestatus[n] == 'non-crispr' ]
    #color and set edge widths
    edgelist = [(u,v) for u,v in net.edges()]
    edge_color = [net[u][v]['weight'] for (u,v) in edgelist]
    edge_width = [net[u][v]['sem_weight']*width_scaling for (u,v) in edgelist]
    #set node positions
    if layout == nx.shell_layout:
        pos = layout(net,[cnodes,ncnodes])
    else:
        pos = layout(net)
    #Draw graph
    cnd = nx.draw_networkx_nodes(net,pos,ax=ax,nodelist=cnodes,node_color=ccol,alpha=alpha)
    ncnd = nx.draw_networkx_nodes(net,pos,ax=ax,nodelist=ncnodes,node_color=nccol,alpha=alpha)
    lbls = nx.draw_networkx_labels(net,pos,nodeidx,ax=ax)
    edges = nx.draw_networkx_edges(net,pos,ax=ax,width=edge_width,edgelist=edgelist,edge_color=edge_color,edge_cmap=edge_cmap)
    #create legends
    pltfnc = partial(plt.legend,fancybox=True,loc='upper right',prop={'size':10})
    fig.colorbar(edges,ax=ax)
    #width legend
    widthlegend = []
    widthticks = np.linspace(min(edge_width),max(edge_width),legend_entries)
    for width in widthticks:
        label = np.around(width/width_scaling,legend_entries)
        widthlegend.append(mpl.lines.Line2D([], [], color='black',linewidth=width,label=label))
    bbox2 = (left_legend_offset,1)
    width_legend = pltfnc(handles=widthlegend,bbox_to_anchor=bbox2,title="$\\bf{Std. Error}$")
    ax = plt.gca().add_artist(width_legend)
    #status legend
    statuslegend = [mpl.patches.Patch(color=ccol,label='CRISPR'),
            mpl.patches.Patch(color=nccol,label='Non-CRISPR')]
    height1 = 1-(len(widthlegend)+1)*entry_width
    bbox1 = (left_legend_offset,height1)
    #status_legend = plt.legend(handles=statuslegend,loc='upper right',bbox_to_anchor=bbox1,title='Node Status',fancybox=True,title_fontsize=16)
    status_legend = pltfnc(handles=statuslegend,bbox_to_anchor=bbox1,title="$\\bf{Node Status}$")
    ax = plt.gca().add_artist(status_legend)
    #accession legend
    acclegends = []
    for a,i in nodeidx.items():
        label = '{} = {}'.format(i,a)
        acclegends.append(mpl.lines.Line2D([],[],linewidth=0,label=label))
    height2 = height1-(len(statuslegend)+1)*entry_width
    bbox3 = (left_legend_offset,height2)
    pltfnc(handles=acclegends,bbox_to_anchor=bbox3,title="$\\bf{Accession Numbers}$")
#    handles = statuslegend + widthlegend + acclegends
#    plt.legend(handles=handles,loc='upper right'
    if type(save) != bool:
        fig.savefig(save,dpi=150,format='png',frameon=False)
    plt.show()

def loadReports():
    jrdir = '/home/sidreed/thesis_SidReed/json_reports/'
    mkdir = '/home/sidreed/thesis_SidReed/marko_reports/'
    markos = [json.load(open(os.path.join(mkdir,path))) for path in os.listdir(mkdir)]
    mdf = pd.DataFrame(markos).set_index('genus')
    stats = [json.load(open(os.path.join(jrdir,path))) for path in os.listdir(jrdir)]
    stats = {s['genus']:s for s in stats}
    commongens = set(mdf.index).intersection(set(list(stats.keys())))
    mdf = mdf.loc[list(commongens),:]
    return mdf,stats

def finaldf(mdf,stats):
    mdf['a'] = abs(mdf['crispr_indel_rate']-mdf['non-crispr_indel_rate'])
    mdf['b'] = [np.mean(stats[i]['modularity']) for i in mdf.index]
    mdf['c'] = [np.mean(stats[i]['assortativity']) for i in mdf.index]
    mdf['d'] = [stats[i]['degree']['crispr_mean'] for i in mdf.index]
    mdf['e'] = [stats[i]['degree']['non-crispr_means'] for i in mdf.index]
    mdf['f'] = [stats[i]['clustering']['crispr_mean'] for i in mdf.index]
    mdf['g'] = [stats[i]['clustering']['non-crispr_means'] for i in mdf.index]
    mdf.columns = ['c_indel','nc_indel','c_sem_indel',
                   'nc_sem_indel','cnc_ratediff','net_mean_mod',
                   'net_mean_assort','c_mean_deg','nc_mean_deg',
                   'c_mean_clust','nc_mean_clust'
                  ]
    return mdf




if __name__ == '__main__':
    print('nothing')
