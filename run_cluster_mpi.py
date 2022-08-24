import glob

from model.cluster import Reducer
import argparse
import pandas as pd
import numpy as np

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

    fps_paths = glob.glob(r'out/paths/*.fps')

    print(fps_paths)

    reducer = Reducer(
        n_clusters=3000,
        batch_size=10000,
        max_iter=1000,
        init_size=3000
    )








