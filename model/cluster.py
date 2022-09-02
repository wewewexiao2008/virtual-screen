import glob
from memory_profiler import profile
import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sklearn.cluster import MiniBatchKMeans
from loguru import logger
from rdkit.Chem import Draw
from scripts.rdkit2pdbqt import MolFromPDBQTBlock
from data.tools.utils import timing, split_n
import datetime
import os
import socket

now = datetime.datetime.now
FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def show_raw_mols(df, mol_dir, layer1, layer2, out_dir):
    mols = []
    ids = []
    for mol_id in df[(df['layer1'] == layer1) & (df['layer2'] == layer2)]:
        _ls = glob.glob(r'{}/**/{}.pdbqt'.format(mol_dir, mol_id), recursive=True)
        mol_path = _ls[0]
        with open(mol_path, 'r') as mol_f:
            mol = MolFromPDBQTBlock(mol_f.read())
        mols.append(mol)
        ids.append(mol_id)

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=4,
        subImgSize=(200, 200),
        legends=[i for i in ids]
    )

    img.save('{}/{}_{}.png'.format(out_dir, layer1, layer2))


def read_base64(df):
    res = []
    for i in df['base64']:
        _fp = ExplicitBitVect(0)
        ExplicitBitVect.FromBase64(_fp, i)
        arr = np.zeros((1,), dtype=np.bool_)
        DataStructs.ConvertToNumpyArray(_fp, arr)
        res.append(arr)

    res = np.array(res)
    res.astype(np.bool_)
    return res


class Reducer:
    """
    Reducer:

    """

    def __init__(self,
                 n_clusters: int = 3000,
                 batch_size: int = 10000,
                 max_iter: int = 1000,
                 init_size: int = 10000,
                 layer: int = 0
                 ):

        self.n_clusters = n_clusters
        self.batch_size = batch_size
        self.max_iter = max_iter
        self.init_size = init_size
        self.layer = layer

    def mb_kmeans(self, X, verbose):

        if verbose:
            logger.info("shape: {}, itemsize: {}".format(X.shape, X.itemsize))

        clustering = MiniBatchKMeans(n_clusters=self.n_clusters,
                                     batch_size=self.batch_size,
                                     max_iter=self.max_iter,
                                     init_size=self.init_size,
                                     init='k-means++', verbose=0, compute_labels=True,
                                     random_state=None, tol=0.0, max_no_improvement=10,
                                     n_init=3, reassignment_ratio=0.01)

        n_samples, n_features = X.shape
        n_batches = int(np.ceil(float(n_samples) / self.batch_size))

        if verbose:
            for i, x_batch in enumerate(split_n(X, n_batches)):
                with timing("partial fitting...{}/{}".format(i, n_batches)):
                    clustering = clustering.partial_fit(x_batch)
        else:
            for x_batch in split_n(X, n_batches):
                clustering = clustering.partial_fit(x_batch)

        if verbose:
            with timing("predicting..."):
                y = clustering.predict(X)
        else:
            y = clustering.predict(X)

        return y, clustering.inertia_

    def run_with_fps_mpi(self, fps_path, out_path, col, verbose):
        if verbose:
            with timing("reading csv"):
                df = pd.read_csv(fps_path, delimiter='\t')
        else:
            df = pd.read_csv(fps_path, delimiter='\t')

        return self.run_with_df_mpi(df, out_path, col, verbose)

    def run_with_df_mpi(self, df, out_path, col, verbose):
        if verbose:
            with timing("reading base64"):
                X = read_base64(df)
        else:
            X = read_base64(df)

        y, inertia = self.mb_kmeans(X, verbose)
        df[col] = y
        joblib.dump(y, out_path)
        return df, inertia

    def run_with_fps(self, fps_path, file_path, out_path, verbose=False):
        """main procedure"""

        with open(file_path + ".log", 'a+') as f:
            f.write('job started on {0}\n'.format(socket.gethostname()))
            f.write('new task for fps_path=' + str(fps_path) + '\n')

        df = pd.read_csv(fps_path, delimiter='\t')
        # with timing("loading csv"):
        if verbose:
            with open(file_path + ".log", 'a+') as f:
                f.write("cluster data amount: {}\n".format(len(df)))

        X = []
        if verbose:
            t0 = now()
            with open(file_path + ".log", 'a+') as f:
                f.write('Start reading as numpy at {0}\n'.format(t0.isoformat()))

        for i in df['base64']:
            _fp = ExplicitBitVect(0)
            ExplicitBitVect.FromBase64(_fp, i)
            arr = np.zeros((1,), dtype=np.bool_)
            DataStructs.ConvertToNumpyArray(_fp, arr)
            X.append(arr)

        if verbose:
            t0 = now()
            with open(file_path + ".log", 'a+') as f:
                f.write('Start fit and predict at {0}\n'.format(t0.isoformat()))

        y, inertia = self.mb_kmeans(X)

        if verbose:
            t1 = now()
            h = (t1 - t0).total_seconds() // 3600
            m = (t1 - t0).total_seconds() // 60 - h * 60
            s = (t1 - t0).total_seconds() - m * 60 - h * 60
            with open(file_path + ".log", 'a+') as f:
                f.write('Finished at {0} after '
                        '{1}h {2}min {3:0.2f}s\n'.format(t1.isoformat(), h, m, s))

        joblib.dump(y, out_path + ".result")

        return fps_path, inertia
