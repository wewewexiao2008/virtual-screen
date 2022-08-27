import glob

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
    parser.add_argument('-n', '--n_cpus', type=int, default=576)

    args = parser.parse_args()

    fps_dir = args.fps_dir
    out_dir = args.out_dir
    n_cpus = args.n_cpus

    reducer = Reducer(
        n_clusters=2000,
        batch_size=1000,
        max_iter=1000,
        init_size=2000
    )

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    root = 100

    if comm_rank == root:
        fps_paths = glob.glob(r'out/fingerprint_mpi_new/*.fps', recursive=True)

        send_buf = fps_paths
    else:
        send_buf = None

    local_data = comm.scatter(send_buf, root=root)

    reducer.run(local_data)





