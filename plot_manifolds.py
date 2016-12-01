#!/usr/bin/python

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, SpectralEmbedding
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

# Read in allelic dosages
dosage_df = pd.read_csv('pruned.raw', sep=' ')

# Import population information
pop_df = pd.read_csv('strs_samples.20130502.ALL.panel', sep='\t', 
    usecols=['sample','pop','super_pop'])
pop_df.columns = ['FID','pop','super_pop']
superpopset_list = list(set(pop_df['super_pop']))

# Merge population assignment info with allelic dosages
df = pd.merge(pop_df,dosage_df, on='FID')