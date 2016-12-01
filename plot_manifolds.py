#!/usr/bin/python

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, SpectralEmbedding
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

from keras.layers import Input, Dense, Lambda
from keras.models import Model
from keras import backend as K
from keras import objectives


# Plots 2D manifold of the transformed input data
def plot_manifold(X_c, superpop_list, superpopset_list, title):
	# Generate separate colors for each superpopulation
	colors = iter(cm.rainbow(np.linspace(0, 1, len(superpopset_list))))
    # Plot samples belonging to each superpopulation with separate colors
	for i in range(len(superpop_list)):
		plt.scatter(X_c[superpop_list[i][0],0], X_c[superpop_list[i][0],1],
			color=next(colors),label=superpopset_list[i])

	plt.xlabel('C1')
	plt.ylabel('C2')
	plt.legend(loc='lower right')
	plt.title(title)
	plt.show()


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

# Principal components
pca = PCA(n_components=2)
X_pc = pca.fit(X).transform(X)
plot_manifold(X_pc, superpop_list, superpopset_list, 'Principal components')

# MDS manifold
mds = MDS(n_components=2)
X_mds = mds.fit_transform(X)
plot_manifold(X_mds, superpop_list, superpopset_list, 'Multidimensional scaling')

# t-SNE
tsne = TSNE(n_components=2)
X_tsne = tsne.fit_transform(X)
plot_manifold(X_tsne, superpop_list, superpopset_list, 't-distributed Stochastic Neighbor Embedding')

# Spectral Embedding
tsne = TSNE(n_components=2)
X_tsne = tsne.fit_transform(X)
plot_manifold(X_tsne, superpop_list, superpopset_list, 'Spectral Embedding')