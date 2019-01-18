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
import plotly.offline as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#drawing
def makecolors(cmap):
    scheme = cm = plt.get_cmap(cmap)
    cNorm = colors.Normalize(vmin=0,vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=scheme)
    return scalarMap.to_rgba

def makeEdge(edge,pos,colormap,maxthick=20):
    x0,y0 = pos[edge[0]]
    x1,y1 = pos[edge[1]]
    color = colormap(edge[2]['weight'])
    newedge = dict(
        type='scatter',
        x = tuple([x0,x1,None]),
        y = tuple([y0,y1,None]),
        hoverinfo='text',
        text = tuple([edge[2]['weight']]),
        mode='lines',
        line=dict(
            width=maxthick*edge[2]['weight'],
            color='rgb{}'.format(color)
            )
        )
    return newedge

def makeNode(node,pos):
    x,y = pos[node[0]]
    color=lambda node:'red' if node[1]['isCRISPR'] else 'blue'
    newnode = dict(type='scatter',
                    x=tuple([x]),
                    y=tuple([y]),
                    text=node[0],
                    mode='markers',
                    hoverinfo='text',
                    marker=dict(
                        color=color(node),
                        size=30,
                   )
              )
    return newnode

def drawNetwork(network,title,colormap=cmx.viridis):
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
                        title=title,
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
    return fig

