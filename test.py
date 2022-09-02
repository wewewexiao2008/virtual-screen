import joblib
import numpy as np
import pandas as pd
from loguru import logger
from memory_profiler import profile
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect
from sklearn.cluster import MiniBatchKMeans
from tqdm import tqdm

from dask.distributed import Client
import joblib

from data.tools.utils import timing

def mb_kmeans(df, verbose=False, save=False):
    fps = []
    for i in tqdm(df['base64']):
        _fp = ExplicitBitVect(0)
        ExplicitBitVect.FromBase64(_fp, i)
        # arr = np.zeros((1,), dtype=np.bool_)
        arr = np.zeros((1,), dtype=np.bool_)
        DataStructs.ConvertToNumpyArray(_fp, arr)
        fps.append(arr)

    if verbose:
        logger.info("fps loaded")
    with joblib.parallel_backend('dask'):
        clustering = MiniBatchKMeans(n_clusters=5,
                                    batch_size=100,
                                    max_iter=100,
                                    init_size=200,
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

    df['layer_{}'.format(1)] = y
    return df


def run_with_fps(fps_path, verbose=False):
    """main procedure"""
    with timing("Reading csv"):
        df = pd.read_csv(fps_path, delimiter='\t')
    if verbose:
        logger.info("cluster data amount: {}".format(len(df)))
    return mb_kmeans(df, verbose=verbose)


if __name__ == "__main__":
    # client = Client(processes=False)  # create local cluster
    fps = 'out/fingerprint_mpi_new/0-cu2.fps'
    # fps = 'packed.fps'

    df = run_with_fps(fps, verbose=True)

    # df.to_csv('test_result_layer1.fps')
