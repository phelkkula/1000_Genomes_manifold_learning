# Manifold learning of 1000 Genomes samples
#### Manifolds using different dimensionality reduction methods 

This repository provides examples for learning manifolds of the 1000 Genomes genotypes using various dimensionality reduction techniques. Plots of the two-dimensional manifolds are in the _plots_ 
folder. At the moment, this repository has examples using the following methods:

1.  Principal Components

2.  Multidimnesional scaling

3.  Spectral Embedding

4.  t-distributed Stochastic Neighbor Embedding

5.  Variational Autoencoder (implementation with Keras)

The manifolds are computed from QC:ed and LD-pruned chr20 variants (see variant_pruning.sh for details). The vcf can be downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
