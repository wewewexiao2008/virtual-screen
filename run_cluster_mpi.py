import os
os.environ['OPENBLAS_NUM_THREADS'] = '12'

import glob
import sys

from loguru import logger

from model.cluster import Reducer
import argparse
import pandas as pd
import numpy as np
import mpi4py.MPI as MPI



def main():

    description = 'Multi layer K-means'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-f', '--fps_dir', type=str, required=True)
    parser.add_argument('-o', '--out_dir', type=str, required=True)

    args = parser.parse_args()

    fps_dir = args.fps_dir
    out_dir = args.out_dir

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    nc_layer1 = comm_size*2
    nc_layer2 = comm_size*4

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

    root = 100
    verbose = True if comm_rank == root else False
    if comm_rank == root:
        log_dir = os.path.join(os.path.curdir, 'log')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'debug_{time}.log')
        logger.add(log_file)

        fps_paths = glob.glob(r'{}/*.fps'.format(fps_dir), recursive=True)
        send_buf = fps_paths
    else:
        send_buf = None

    # path list
    local_data = comm.scatter(send_buf, root=root)

    """Layer1"""
    sys.stdout.write("Running Layer1...\n")
    df_raw = pd.read_csv(local_data, delimiter='\t')
    if verbose:
        logger.info("cluster data amount: {}".format(len(df_raw)))
    df = l1_reducer.run_with_fps(df_raw, verbose)

    try:
        df.to_csv('{}/layer1/ckpt_{}.tsv'.format(out_dir, comm_rank), sep='\t')
    except Exception as e:
        os.makedirs('{}/layer1'.format(out_dir))
        df.to_csv('{}/layer1/ckpt_{}.tsv'.format(out_dir, comm_rank), sep='\t')

    """Layer2"""
    sys.stdout.write("Running Layer2...\n")
    l2_dfs = []
    for i in range(nc_layer1):
        _df = df[df['layer_1'] == i]
        _df = l2_reducer.run_with_df(_df, verbose)
        l2_dfs.append(_df)

    df_final = pd.concat(l2_dfs)
    try:
        df_final[['id', 'layer_1', 'layer_2']].to_csv('{}/layer2/ckpt_{}.tsv'.format(out_dir, comm_rank), sep='\t')
    except Exception as e:
        os.makedirs('{}/layer2'.format(out_dir))
        df_final[['id', 'layer_1', 'layer_2']].to_csv('{}/layer2/ckpt_{}.tsv'.format(out_dir, comm_rank), sep='\t')




if __name__ == "__main__":
    main()


