import os

import joblib

os.environ['OPENBLAS_NUM_THREADS'] = '60'

import glob
import sys

from loguru import logger

from model.cluster import Reducer
import argparse
import pandas as pd
import numpy as np
import mpi4py.MPI as MPI
from data.tools.utils import split_n, timing

FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    description = 'Multi layer K-means'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--fps_dir', type=str, required=True)
    parser.add_argument('-o', '--out_dir', type=str, required=True)
    parser.add_argument('-l', '--log_dir', type=str)

    args = parser.parse_args()

    fps_dir = args.fps_dir
    out_dir = args.out_dir
    log_dir = args.log_dir

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    nc_layer1 = 10000
    nc_layer2 = 200

    l1_reducer = Reducer(
        n_clusters=nc_layer1,
        batch_size=1024,
        max_iter=1000,
        init_size=nc_layer1,
        layer=1
    )
    l2_reducer = Reducer(
        n_clusters=nc_layer2,
        batch_size=1000,
        max_iter=1000,
        init_size=nc_layer2,
        layer=2
    )

    root = 0
    verbose = True if comm_rank == root else False
    if comm_rank == root:
        logger.add(os.path.join(log_dir, 'debug.log'))
        logger.info("comm_size = {}".format(comm_size))
        fps_paths = glob.glob(r'{}/*.fps'.format(fps_dir), recursive=True)
        send_buf = [i for i in split_n(fps_paths, comm_size)]
    else:
        send_buf = None

    # path list
    local_data = comm.scatter(send_buf, root=root)

    for i, fps_path in enumerate(local_data):
        """Layer1"""
        basename = os.path.basename(fps_path)
        out_path = os.path.join(out_dir, basename + '.layer1')

        if verbose:
            with timing("running {} layer1".format(basename)):
                df, inertia = l1_reducer.run_with_fps_mpi(fps_path, out_path, 'layer1', verbose)
        else:
            df, inertia = l1_reducer.run_with_fps_mpi(fps_path, out_path, 'layer1', verbose)

        """Layer2"""
        l2_dfs = []
        stat_ls = []
        # layer 1 inertia
        for j in range(nc_layer1):
            if verbose and j % 100 == 0:
                logger.info("Layer2: {}/{} done".format(j, nc_layer1))

            _df = df[df['layer1'] == j]
            out_path = os.path.join(out_dir, basename + '.' + str(j) + '.layer2')
            _df, inertia = l2_reducer.run_with_df_mpi(_df, out_path, 'layer2', verbose)
            l2_dfs.append(_df)
            stat_ls.append(inertia)

        if comm_rank == root:
            logger.info("output stat to {}".format(basename + '.stat'))

        stat_df = pd.DataFrame(stat_ls, columns=['inertia'])
        stat_df.to_csv(os.path.join(out_dir, basename + '.stat'))

        df_final = pd.concat(l2_dfs)

        df_final[['id', 'layer1', 'layer2']].to_csv('{}/final_{}.tsv'.format(out_dir, basename), sep='\t')


if __name__ == "__main__":
    main()
