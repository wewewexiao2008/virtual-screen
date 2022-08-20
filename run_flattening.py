import os
from loguru import logger
from data.tools import utils, pipeline
import sys
import glob
import shutil
import multiprocessing


def split_n(origin_list, n):
    l = len(origin_list)
    import math
    res = [origin_list[math.floor(i / n * l):math.floor((i + 1) / n * l)] for i in range(n)]
    for i in res:
        yield i


def flatten(path, tmp_dir):
    # filename = os.path.basename(path)
    logger.info(path)
    shutil.move(path, tmp_dir)


def main():
    import argparse

    description = 'flatten pdbqt files.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--tmp_dir', type=str, required=True, help='input directory')

    args = parser.parse_args()

    tmp_dir = args.tmp_dir  # './tmp/tmp6sj_2ixi'

    with utils.timing("counting pdbqt files"):
        sys.stdout.write("counting pdbqt files...\n")
        paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)

    with utils.timing("flattening"):

        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        logger.info(multiprocessing.cpu_count())
        for p in paths:
            pool.apply_async(flatten, args=(p, tmp_dir))
        pool.close()
        pool.join()

    shutil.rmtree()


if __name__ == "__main__":
    main()
