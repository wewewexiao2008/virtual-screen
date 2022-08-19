import multiprocessing
import os
from loguru import logger
from data.tools import utils, pipeline
from model.cluster import Reducer


def main():
    import argparse

    description = 'fingerprinting from pdbqt files.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-o', '--out_path', type=str, metavar='<out_path>', required=True,
                        help='output directory')
    parser.add_argument('-n', '--n_cpu', type=int,
                        help='assign the number of processes.')

    args = parser.parse_args()

    out_path = args.out_path
    n_cpu = args.n_cpu

    log_dir = os.path.join(os.path.curdir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file)

    if not os.path.exists('./out/'):
        os.makedirs('./out/')

    data_pipeline = pipeline.DataPipeline(
        data_dir='.',
        n_cpu=n_cpu,
        fp_type='morgan'
    )

    tmp_dir = './tmp/tmp6sj_2ixi'
    fps_path = os.path.join('./out/', out_path)

    with utils.timing("calculating fingerprint"):
        data_pipeline.fingerprint(tmp_dir, fps_path)


if __name__ == "__main__":
    main()
