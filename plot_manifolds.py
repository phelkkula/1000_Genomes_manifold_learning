#!/usr/bin/python

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, SpectralEmbedding
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

# Read in allelic dosages
dosage_df = pd.read_csv('pruned.raw', sep=' ')