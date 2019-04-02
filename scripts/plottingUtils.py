import os
import json
import numpy as np
import pandas as pd
import networkx as nx
import scipy.stats as sst
from functools import partial
from functools import reduce
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.inset_locator as ins

#Utils
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

def cnodescount(genus,andf):
    gdf = andf[andf['Genus'] == genus].fillna(0)
    gdf = gdf[gdf['nucPath'] != 0]
    t = gdf.shape[0]
    c = gdf[gdf['isCRISPR'] == True].shape[0]
    return c/t,t

def numCrisprNodes(genus):
    path = '/home/sidreed/thesis_SidReed/plotnets/{}_plot_network.json'.format(genus)
    df = pd.DataFrame(json.load(open(path)))
    df = df[df['source'] != df['sink']].drop('target',axis=1)
    status = getNodeStatus(df)
    return sum([1 for k,v in status.items() if v == 'crispr']),len(list(status.keys()))

def genlist(stats,col,thresh):
    #ratios = [k for k,v in stats.items() if np.mean(v[col])/sst.sem(v[col]) < 100]
    ratios = [(k,abs(np.mean(v[col])/sst.sem(v[col])))  for k,v in stats.items()]
    ratios = [r[0] for r in ratios if r[1] < thresh]
    return ratios


#Data IO
def loadnetwork(path,dfonly=False,returnname=False):
    netjson = json.load(open(path))
    df = pd.DataFrame(netjson)
    df = df[df['source'] != df['sink']].drop('target',axis=1)
    if not dfonly:
        net = nx.from_pandas_edgelist(df,'source','sink',edge_attr=True)
        if returnname:
            name = os.path.basename(path).split('_')[0]
            return df,net,name
        else:
            return df,net
    else:
        if returnname:
            name = os.path.basename(path).split('_')[0]
            return df,name
        else:
            return df

def loadMarkoOnly(andf):
    mkdir = '/home/sidreed/thesis_SidReed/marko_reports'
    mkdir = '/home/sidreed/thesis_SidReed/all_markos/'
    markos = [json.load(open(os.path.join(mkdir,path))) for path in os.listdir(mkdir)]
    mdf = pd.DataFrame(markos)
    mdf = mdf[pd.notnull(mdf['genus'])]
    ccount = partial(cnodescount,andf=andf)
    counts = mdf['genus'].apply(ccount)
    mdf['a'],mdf['b'] = [x[0] for x in counts],[x[1] for x in counts]
    mdf = mdf.set_index('genus')
    mdf.columns = ['c_indel','c_sem_indel','nc_indel','nc_sem_indel','c_otus','t_otus']
    mdf = mdf[(mdf['c_indel'] != 100) & (mdf['nc_indel'] != 100)]
    return mdf

def loadReports():
    mkdir = '/home/sidreed/thesis_SidReed/marko_reports/'
    markos = [json.load(open(os.path.join(mkdir,path))) for path in os.listdir(mkdir)]
    mdf = pd.DataFrame(markos).set_index('genus')
    jrdir = '/home/sidreed/thesis_SidReed/json_reports/'
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
    c_t_counts = [numCrisprNodes(i) for i in mdf.index]
    mdf['h'] = [x[1] for x in c_t_counts]
    mdf['i'] = [x[0]/x[1] for x in c_t_counts]
    mdf.columns = ['c_indel','c_sem_indel','nc_indel',
                   'nc_sem_indel','cnc_ratediff','net_mean_mod',
                   'net_mean_assort','c_mean_deg','nc_mean_deg',
                   'c_mean_clust','nc_mean_clust','t_otus','c_otus'
                  ]
    return mdf

#Plotting
def plotNetwork(net,df,width_scaling=4000,alpha=0.8,legend_entries=6,legend_decimals=6,ccol='blue',nccol='red',edge_cmap=plt.cm.viridis_r,layout=nx.spring_layout,seed=9,entry_width=0.03*0.63,left_legend_offset=1.05,dpi=50,name='',save=False):
    #set up plot
    np.random.seed(seed)
    fig,ax = plt.subplots(1,1,figsize=(20,15))
    ax.set_axis_off()
    ax.set_title(name)
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
    #savefig
    if type(save) != bool:
        fig.savefig(save,dpi=dpi,format='png',frameon=False)
    plt.show()

def multiBarPlot(df,cols,ylabel,xlabel='Genera',width=1,dpi=50,
                labels=['CRISPR','Non-CRISPR'],file=False):
    sdf = df.sort_values(by=cols)
    fig, ax = plt.subplots(figsize=(20,10))
    pos = np.arange(0,len(sdf[cols[0]])*(len(cols)+1),len(cols)+1)
    pal = sns.color_palette('coolwarm')
    colors = [pal[0],pal[-1]]
    for i,col in enumerate(cols):
        ax.bar(pos+width*i,sdf[cols[i]],width,
                color=colors[i],label=labels[i])
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_xticks(pos+0.5*width)
    ax.set_xticklabels(sdf.index,rotation=65,ha='right')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend()
    if type(file) != bool:
        fig.savefig(file,dpi=dpi,format='png',frameon=False)
    plt.show()

def cVsncRate(nohdf,size=(10,5),dpi=50,file=False):
    fig, ax = plt.subplots(figsize=size)
    #data
    x,y = 'nc_indel', 'c_indel'
    xerr = nohdf['nc_sem_indel']
    yerr = nohdf['c_sem_indel']
    #main fig
    sns.scatterplot(x=x,y=y,data=nohdf,ax=ax)
    #ax.errorbar(nohdf[x], nohdf[y], xerr=xerr,yerr=yerr, fmt='o',ecolor='lightgray')
    ax.set_xlabel('Non-CRISPR Gene Indel Rate')
    ax.set_ylabel('CRISPR Gene Indel Rate')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #inlet fig
    axins = ins.zoomed_inset_axes(ax,4.5,loc=1)
    sns.scatterplot(x=x,y=y,data=nohdf,ax=axins)
    #axins.errorbar(nohdf[x], nohdf[y], xerr=xerr,yerr=yerr, fmt='o',ecolor='lightgray')
    axins.set_xlim(0,10.5)
    axins.set_ylim(-2,15)
    axins.xaxis.set_ticks_position('none')
    axins.yaxis.set_ticks_position('none')
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.set_xlabel('')
    axins.set_ylabel('')
    #axins.spines['right'].set_visible(False)
    #axins.spines['top'].set_visible(False)
    ins.mark_inset(ax,axins,loc1=2,loc2=4,fc='none',ec='0.5')
    #wilcoxon annotate
    wilx = sst.wilcoxon(nohdf[x],nohdf[y],zero_method='pratt')
    text = 'Wilcoxon Rank: {}\nP-Value: {}'.format(wilx.statistic,
                                                np.round(wilx.pvalue,10))
    axins.annotate(text,xy=(0.05,0.80),xycoords='axes fraction')
    if type(file) != bool:
        fig.savefig(file,dpi=dpi,format='png',frameon=False)
    plt.show()

def cVsncClust(nohdf,size=(10,5),dpi=50,file=False):
    fig, ax = plt.subplots(figsize=size)
    cols = sns.color_palette('coolwarm')
    blue,red = cols[0],cols[-1]
    #lin reg
    x,y = 'c_mean_clust', 'nc_mean_clust'
    xdata = nohdf[x].astype(float).values
    ydata = nohdf[y].astype(float).values
    m,b,r,p,se = sst.linregress(xdata,ydata)
    y = lambda x:m*x+b
    #main fig
    ax.scatter(xdata,ydata,color=blue)
    ax.plot(xdata,y(xdata),color=red,
       label='$R^2$: {}  P-val: {}'.format(np.round(r**2,5),np.round(p,20)))
    ax.set_xlabel('Non-CRISPR Mean Clustering Coefficient')
    ax.set_ylabel('CRISPR Mean Clustering Coefficient')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend()
    #inlet fig
    axins = ins.zoomed_inset_axes(ax,1.55,loc=4)
    axins.scatter(xdata,ydata,color=blue)
    axins.plot(xdata,y(xdata),color=red)
    axins.set_xlim(0.03,0.16)
    axins.set_ylim(0.02,0.17)
    axins.xaxis.set_ticks_position('none')
    axins.yaxis.set_ticks_position('none')
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.set_xlabel('')
    axins.set_ylabel('')
    ins.mark_inset(ax,axins,loc1=2,loc2=4,fc='none',ec='0.5')
    #wilcoxon annotate
    wilx = sst.wilcoxon(xdata,ydata,zero_method='pratt')
    text = 'Wilcoxon Rank: {}\nP-Value: {}'.format(wilx.statistic,
                                                np.round(wilx.pvalue,5))
    axins.annotate(text,xy=(0.43,0.05),xycoords='axes fraction')
    if type(file) != bool:
        fig.savefig(file,dpi=dpi,format='png',frameon=False)
    plt.show()

def rateVsCFrac(nohdf,size=(10,5),dec=5,dpi=50,file=False):
    #Data
    x = nohdf['c_otus']
    y1 = nohdf['c_indel']
    y2 = nohdf['nc_indel']
    #base plotting
    fig,ax = plt.subplots(figsize=size)
    cols = sns.color_palette('coolwarm')
    ax.scatter(x=x,y=y1,color=cols[0],marker='+',label='CRISPR')
    ax.scatter(x=x,y=y2,color=cols[-1],marker='x',label='Non-CRISPR')
    #CRISPR
    m,b,r,p,se = sst.linregress(x,y1)
    y = lambda x:m*x+b
    ax.plot(x,y(x),color=cols[0],label='$R^2$: {} P-value: {}'.format(np.round(r**2,dec),np.round(p,dec)))
    #Non-CRISPR
    m,b,r,p,se = sst.linregress(x,y2)
    y = lambda x:m*x+b
    ax.plot(x,y(x),color=cols[-1],label='$R^2$: {} P-value: {}'.format(np.round(r**2,dec),np.round(p,dec)))
    #plot formatting
    plt.legend()
    ax.set_xlabel('Fraction of OTUs With A CRISPR System')
    ax.set_ylabel('Gene Indel Rate of OTUs')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if type(file) != bool:
        fig.savefig(file,dpi=dpi,format='png',frameon=False)
    plt.show()

def violinPlot(stats,col='modularity',thresh=np.inf,
                size=(20,10),dpi=50,limit=-1,file=False):
    reps = len(list(stats.values())[0]['modularity'])
    s1 = genlist(stats,col,thresh)
    s2 = [[s]*reps for s in s1]
    s3 = reduce((lambda x,y:x+y),s2)
    v1 = [stats[g][col] for g in s1]
    v2 =  reduce((lambda x,y:x+y),v1)
    df = pd.DataFrame([s3,v2]).transpose()
    df.columns = ['genera',col]
    df[col] = df[col].astype(float)
    fig,ax = plt.subplots(figsize=size)
    sns.violinplot(y=col,x='genera',data=df,ax=ax,
                    inner='quart',scale='count')
    ax.set_xticklabels(s1,rotation=45,ha='right')
    ax.set_xlabel('Genera')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('{}{}'.format(col[0].upper(),col[1:]))
    if type(file) != bool:
        fig.savefig(file,dpi=dpi,format='png',frameon=False)
    plt.show()
    return ax


if __name__ == '__main__':
    print('nothing')

#DEPRECIATED
#def violinModAsstPlot():
#    #make df
#    jrdir = '/home/sidreed/thesis_SidReed/json_reports/'
#    stats = [json.load(open(os.path.join(jrdir,path))) for path in os.listdir(jrdir)]
#    stats = {s['genus']:s for s in stats}
#    reps = len(list(stats.values())[0]['modularity'])
#    genera = list(stats.keys())
#    ngenera = len(genera)
#    statns = ['modularity','assortativity']
#    ls = len(statns)
#    glist = combl([[g]*ls*reps for g in genera])
#    slist = combl([[stat]*ngenera*reps for stat in statns])
#    vlist = combl([stats[g]['modularity']+stats[g]['assortativity'] for g in genera])
#    df = pd.DataFrame([glist,slist,vlist]).transpose()
#    df.columns = ['genera','stat','value']
#    df.value = df.value.astype(float)
#    #make plot
#    fig,ax = plt.subplots(figsize=(15,30))
#    sns.violinplot(x="value",y="genera",hue='stat',
#                    split=True,data=df,ax=ax)
#    return df,ax
