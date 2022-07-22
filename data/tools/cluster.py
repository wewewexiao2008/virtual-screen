from rdkit.DataStructs import ExplicitBitVect
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.Fingerprints.ClusterMols import GetDistanceMatrix
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

class Cluster
def kmeans(self, fps: list):
    kmeans = KMeans(n_clusters=5, n_jobs=-1)
    kmeans.fit(fps)
    # PCA
    pca = PCA(n_components=2)
    decomp = pca.fit_transform(fps)

def maxmin(self, fps: list):
    def distij(i, j, fps):
        return 1 - DataStructs.DiceSimilarity(fps[i], fps[j])
    picker = MaxMinPicker()
    pickIndices = picker.LazyPick(distij, len(fps), 10, seed=23)
    # picks = [mol1[x] for x in pickIndices]
