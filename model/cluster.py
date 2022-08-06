import datetime
import time

import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sklearn.cluster import KMeans, MiniBatchKMeans, DBSCAN, AgglomerativeClustering
from sklearn.decomposition import PCA
from loguru import logger
from sklearn.manifold import TSNE

from data.tools.utils import timing


class Reducer:
    def __init__(self,
                 data_path: str,
                 alg: str,
                 dim_reducer: str):
        self.data_path = data_path
        self.alg = alg
        self.dim_reducer = dim_reducer
        self.data_df = pd.DataFrame()
        self.fps = []

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
            clustering = KMeans(n_clusters=8, init='k-means++', n_init=10,
                                max_iter=300, tol=1e-4, precompute_distances='auto',
                                verbose=0, random_state=None, copy_x=True,
                                n_jobs=None,algorithm='auto')
        elif self.alg == 'mb-kmeans':
            clustering = MiniBatchKMeans(n_clusters=8, init='k-means++', max_iter=100,
                                         batch_size=100, verbose=0, compute_labels=True,
                                         random_state=None, tol=0.0, max_no_improvement=10,
                                         init_size=None, n_init=3, reassignment_ratio=0.01)
        elif self.alg == 'dbscan':
            clustering = DBSCAN(eps=0.5, min_samples=5, metric='euclidean',
                                metric_params=None, algorithm='auto', leaf_size=30, p=None,
                                n_jobs=None)
        elif self.alg == 'ward':
            clustering = AgglomerativeClustering(n_clusters=8, linkage='ward')
        else:
            # error handle
            return

        """step 3: fit the data"""
        with timing("fitting with {}".format(self.alg)):
            tic = time.time()
            clustering.fit(self.fps)
            toc = time.time()

        t_train = toc-tic

        """step 4: save the model"""
        joblib.dump(clustering,
                    'checkpoint_{}_{}.pkl'.format(self.alg, datetime.datetime.now().__str__()))

        """step 5: visualize(optional)"""
        if self.dim_reducer == 'tsne':
            reducer = TSNE()


        """step 6: reduce"""
        # TODO: filter strategy

