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

# Extract sample ID's
sample_list = list(dosage_df.FID)

# Extract only allelic dosages
dosage_df = df.iloc[:,8:-1]

# Drop variants with a missing genotype
dosage_df = dosage_df.dropna(axis=1)
# Input allelic dosages
X = dosage_df.as_matrix()

# Get sample indices for each super population
superpop_list = [np.where(list(df['super_pop'] == superpopset_list[i]))
                for i in range(len(superpopset_list))]