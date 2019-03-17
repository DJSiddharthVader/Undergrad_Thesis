import re
import os
import sys
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms import approximation as aprx
from scipy.stats import mannwhitneyu
import matplotlib.cm as cmx
import plotly.offline as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.colors as colors

crisprdfpath = '/home/sid/thesis_SidReed/data/fasta_crispr_annotation_df.json'

def loaddf(csv,andf):
    df = pd.DataFrame(path,sep='~')

def loadAll(cdir):
    return [loaddf(csv) for csv in os.listdir(cdir)]

def avgEdgeWeight(netlist):
    means = [np.mean(net['weight']) for net in netlist]
