import glob

import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sklearn.cluster import MiniBatchKMeans
from loguru import logger
from rdkit.Chem import Draw
from scripts.rdkit2pdbqt import MolFromPDBQTBlock
from data.tools.utils import timing


def show_raw_mols(df, mol_dir, group_id, out_dir):
    mols = []
    ids = []
    for mol_id in df[df['layer_1'] == group_id]['id'][:]:
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

    img.save('{}/{}'.format(out_dir, group_id))


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

    def mb_kmeans(self, df, verbose=False, save=False):
        fps = []
        for i in df['base64']:
            _fp = ExplicitBitVect(0)
            ExplicitBitVect.FromBase64(_fp, i)
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(_fp, arr)
            fps.append(arr)

        if verbose:
            logger.info("fps loaded")

        clustering = MiniBatchKMeans(n_clusters=self.n_clusters,
                                     batch_size=self.batch_size,
                                     max_iter=self.max_iter,
                                     init_size=self.init_size,
                                     init='k-means++', verbose=0, compute_labels=True,
                                     random_state=None, tol=0.0, max_no_improvement=10,
                                     n_init=3, reassignment_ratio=0.01)

        if verbose:
            with timing("running mb-kmeans"):
                y = clustering.fit_predict(fps)
        else:
            y = clustering.fit_predict(fps)

        if save:
            joblib.dump(clustering, '{}'.format())

        df['layer_{}'.format(self.layer)] = y
        return df

    def run_with_fps(self, fps_path, verbose=False):
        """main procedure"""
        df = pd.read_csv(fps_path, delimiter='\t')
        if verbose:
            logger.info("cluster data amount: {}".format(len(df)))
        return self.mb_kmeans(df, verbose=verbose)

    def run_with_df(self, df, verbose=False):
        if verbose:
            logger.info("cluster data amount: {}".format(len(df)))
        return self.mb_kmeans(df, verbose=verbose)