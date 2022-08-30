import os

os.environ['OPENBLAS_NUM_THREADS'] = '60'
import glob
import sys

from loguru import logger

from model.cluster import Reducer
import argparse
import pandas as pd
import numpy as np
from data.tools.utils import split_n
from dask.distributed import Client
import joblib


def main():
    client = Client(processes=False)  # create local cluster
    description = 'Multi layer K-means'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--fps_dir', type=str, required=True)
    parser.add_argument('-o', '--out_dir', type=str, required=True)

    args = parser.parse_args()

    fps_dir = args.fps_dir
    out_dir = args.out_dir

    nc_layer1 = 1000
    nc_layer2 = 2000

    l1_reducer = Reducer(
        n_clusters=nc_layer1,
        batch_size=1000,
        max_iter=1000,
        init_size=nc_layer1
    )
    l2_reducer = Reducer(
        n_clusters=nc_layer2,
        batch_size=1000,
        max_iter=1000,
        init_size=nc_layer2
    )

    root = 0
    verbose = True

    log_dir = os.path.join(os.path.curdir, 'log')
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file)

    fps_paths = glob.glob(r'{}/*.fps'.format(fps_dir), recursive=True)

    with joblib.parallel_backend('dask'):
        for i, fps_path in enumerate(fps_paths):
            """Layer1"""
            # if comm_rank == root:
            sys.stdout.write("{}:{} running Layer1...\n".format(i, fps_path))
            df = l1_reducer.run_with_fps(fps_path, verbose)
            try:
                df.to_csv('{}/layer1/ckpt_{}.tsv'.format(out_dir, i), sep='\t')
            except Exception as e:
                os.makedirs('{}/layer1'.format(out_dir))
                df.to_csv('{}/layer1/ckpt_{}.tsv'.format(out_dir, i), sep='\t')

            """Layer2"""
            sys.stdout.write("{}:{} running Layer2...\n".format(i, fps_path))
            l2_dfs = []
            for j in range(nc_layer1):
                _df = df[df['layer_1'] == j]
                _df = l2_reducer.run_with_df(_df, verbose)
                l2_dfs.append(_df)

            df_final = pd.concat(l2_dfs)
            try:
                df_final[['id', 'layer_1', 'layer_2']].to_csv('{}/layer2/ckpt_{}.tsv'.format(out_dir, i), sep='\t')

            except Exception as e:
                os.makedirs('{}/layer2'.format(out_dir))
                df_final[['id', 'layer_1', 'layer_2']].to_csv('{}/layer2/ckpt_{}.tsv'.format(out_dir, i), sep='\t')


if __name__ == "__main__":
    main()
