import tarfile
import os
from loguru import logger
from tqdm import tqdm
from openbabel import openbabel
import concurrent.futures
from data.tools import utils, pipeline
# from data.tools import extractor


def main():
    import argparse

    description = 'A simple command-line interface for virtual flow.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--data_dir', type=str, metavar='<data_path>', required=True,
                        help='input directory')
    parser.add_argument('-o', '--out_dir', type=str, metavar='<out_path>', required=True,
                        help='output directory')
    parser.add_argument('-n', '--n_cpu', type=int, default=16,
                        help='assign the number of processes.')
    parser.add_argument('-t', '--type', type=str, default='morgan',
                        help='fingerprint type: morgan or maccs.')
    # group = parser.add_mutually_exclusive_group(required=True)

    args = parser.parse_args()

    """config"""
    data_dir = args.data_dir
    out_dir = args.out_dir
    n_cpu = args.n_cpu
    fp_type = args.type

    data_pipeline = pipeline.DataPipeline(
        data_dir=data_dir,
        out_dir=out_dir,
        n_cpu=n_cpu,
        fp_type=fp_type
    )

    # logging
    log_dir = os.path.join(os.path.curdir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file)

    """main procedure"""
    with utils.tmpdir_manager('.', delete=True) as tmp_dir:
        with utils.timing("extracting and calculating fingerprint"):
            data_pipeline.extract(tmp_dir)
            data_pipeline.fingerprint(tmp_dir)


if __name__ == "__main__":
    main()
