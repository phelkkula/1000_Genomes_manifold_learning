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

'''
Kingma DK, Welling M "Auto-Encoding Variational Bayes" https://arxiv.org/abs/1312.6114
Variational Autoencoder code adapted from:
https://github.com/fchollet/keras/blob/master/examples/variational_autoencoder.py
'''

from keras.layers import Input, Dense, Lambda
from keras.models import Model
from keras import backend as K
from keras import objectives

# Batch size has to be a product of the prime factors of the sample size
batch_size = 313
original_dim = X.shape[1]
latent_dim = 2
intermediate_dim = 3378
nb_epoch = 100

x = Input(batch_shape=(batch_size, original_dim))
h = Dense(intermediate_dim, activation='relu')(x)
z_mean = Dense(latent_dim)(h)
z_log_var = Dense(latent_dim)(h)


def sampling(args):
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=(batch_size, latent_dim), mean=0.)
    return z_mean + K.exp(z_log_var / 2) * epsilon

z = Lambda(sampling, output_shape=(latent_dim,))([z_mean, z_log_var])

# we instantiate these layers separately so as to reuse them later
decoder_h = Dense(intermediate_dim, activation='relu')
decoder_mean = Dense(original_dim, activation='sigmoid')
h_decoded = decoder_h(z)
x_decoded_mean = decoder_mean(h_decoded)


def vae_loss(x, x_decoded_mean):
    xent_loss = original_dim * objectives.binary_crossentropy(x, x_decoded_mean)
    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
    return xent_loss + kl_loss

vae = Model(x, x_decoded_mean)
vae.compile(optimizer='rmsprop', loss=vae_loss)


x_train = X
x_test = x_train
y_test = df.super_pop

vae.fit(x_train, x_train,
        shuffle=True,
        nb_epoch=nb_epoch,
        batch_size=batch_size,
        validation_data=(x_test, x_test))

# build a model to project inputs on the latent space
encoder = Model(x, z_mean)

# Transform input allelic dosages to the learned latent space
x_test_encoded = encoder.predict(x_test, batch_size=batch_size)

# Plot the two-dimensional VAE manifold
plot_manifold(x_test_encoded, superpop_list, superpopset_list, 'Variational Autoencoder')