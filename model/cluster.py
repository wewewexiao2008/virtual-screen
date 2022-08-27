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
    """
    Reducer:

    """
    def __init__(self,
                 n_clusters: int = 3000,
                 batch_size: int = 10000,
                 max_iter: int = 1000,
                 init_size: int = 10000,
                 ):

        self.n_clusters = n_clusters
        self.batch_size = batch_size
        self.max_iter = max_iter
        self.init_size = init_size

    def run(self, fps_path):
        """main procedure"""

        """step 1: read data"""
        fps = []
        df = pd.read_csv(fps_path, delimiter='\t')
        logger.info("cluster data amount: {}".format(len(df)))
        for i in df['base64']:
            _fp = ExplicitBitVect(0)
            ExplicitBitVect.FromBase64(_fp, i)
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(_fp, arr)
            fps.append(arr)

        """step 2: configure"""
        clustering = MiniBatchKMeans(n_clusters=self.n_clusters,
                                     batch_size=self.batch_size,
                                     max_iter=self.max_iter,
                                     init_size=self.init_size,
                                     init='k-means++', verbose=0, compute_labels=True,
                                     random_state=None, tol=0.0, max_no_improvement=10,
                                     n_init=3, reassignment_ratio=0.01)

        """step 3: fit the data"""
        with timing("fitting and predicting with mini-batch k-means"):
            tic = time.time()
            y = clustering.fit_predict(fps)
            toc = time.time()

        # label
        df_new = pd.DataFrame()
        df_new['id'] = df['id']
        df_new['gp1'] = y
        df_new.to_csv(fps_path)
        t_train = toc-tic




