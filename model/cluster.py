import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA
from loguru import logger
import config

from data.tools.utils import timing


class Reducer:
    def __init__(self,
                 data_path: str,
                 alg: str,):
        self.data_path = data_path
        self.alg = alg
        self.data_df = pd.DataFrame()
        self.fps = []

    def kmeans(self, fps: list):
        kmeans = KMeans(n_clusters=5, n_jobs=-1)
        kmeans.fit(fps)
        # PCA
        pca = PCA(n_components=2)
        decomp = pca.fit_transform(fps)

    def run(self):
        """main procedure"""

        """step 1: read data"""
        df = pd.read_csv(self.data_path, delimiter='\t')
        logger.info("cluster data amount: {}".format(len(df)))
        for i in df['base64']:
            _fp = ExplicitBitVect(0)
            ExplicitBitVect.FromBase64(_fp, i)
            self.fps.append(_fp)

        """step 2: configure"""
        if self.alg == 'kmeans':
            # TODO Config Dict module
            clustering = KMeans(
                n_clusters=8, init='k-means++', n_init=10,
                max_iter=300, tol=1e-4, precompute_distances='auto',
                verbose=0, random_state=None, copy_x=True,
                n_jobs=None,algorithm='auto')
        elif self.alg == 'mb-kmeans':
            clustering = MiniBatchKMeans(
                n_clusters=8, init='k-means++', max_iter=100,
                batch_size=100, verbose=0, compute_labels=True,
                random_state=None, tol=0.0, max_no_improvement=10,
                init_size=None, n_init=3, reassignment_ratio=0.01)
        else:
            # error handle
            return

        """step 3: fit the data"""
        with timing("fitting with {}".format(self.alg)):
            clustering.fit(self.fps)

        """step 4: save the model"""

        """step 5: visualize(optional)"""

