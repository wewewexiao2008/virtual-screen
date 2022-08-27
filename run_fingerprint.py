import os
from loguru import logger
from data.tools import utils, pipeline
import mpi4py.MPI as MPI
import sys
import glob
import shutil
import multiprocessing
import joblib
from tqdm import tqdm


def split_n_old(origin_list, n):
    res = [[] for i in range(n)]
    for i, k in enumerate(origin_list):
        res[i % n].append(k)
    for i in res:
        yield i


def split_n(origin_list, n):
    l = len(origin_list)
    import math
    res = [origin_list[math.floor(i / n * l):math.floor((i + 1) / n * l)] for i in range(n)]
    for i in res:
        yield i


def flatten(path, tmp_dir):
    shutil.move(path, tmp_dir)


def saving(blk_id, blk):
    with open("paths_blk{}.txt".format(blk_id), 'wb') as wf:
        joblib.dump(blk, wf, protocol=4, compress=3)


def main():
    import argparse

    description = 'fingerprinting from pdbqt files.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--tmp_dir', type=str, required=True, help='input directory')
    parser.add_argument('-o', '--out_dir', type=str, help='fps output dir', default='./out/')
    parser.add_argument('-n', '--num_block', type=int, default=10)

    args = parser.parse_args()

    out_dir = args.out_dir
    tmp_dir = args.tmp_dir  # './tmp/tmp6sj_2ixi'
    n_blk = args.num_block

    data_pipeline = pipeline.DataPipeline()

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    root = 100

    saving_flag = False

    # manager process
    if comm_rank == root:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        log_dir = os.path.join(os.path.curdir, 'log')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'debug_{time}.log')
        logger.add(log_file)

        if saving_flag:
            with utils.timing("counting pdbqt files"):
                sys.stdout.write("counting pdbqt files...\n")
                paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)
            blocks = [blk for blk in split_n(paths, n_blk)]
            with utils.timing("saving list"):
                for blk_id, blk in tqdm(enumerate(blocks)):
                    saving(blk_id, blk)
            sys.stdout.write("mol num:{}, blk num:{}\n".format(len(paths), n_blk))
            logger.info("mol num:{}".format(len(paths)))
        # else:

        send_buf = ['./out/paths/paths_blk{}.txt'.format(i) for i in range(n_blk)]
        logger.info("blk num:{}".format(len(send_buf)))

    else:
        send_buf = None

    local_data = comm.scatter(send_buf, root=root)

    fps_path = os.path.join(out_dir, '{}-{}.fps'.format(comm_rank, proc_name))
    with open(fps_path, 'w') as wf:
        wf.writelines("id\tbase64\n")

    sys.stdout.write("process {} of {} on {}, handling {} mols, to {}\n".format(
        comm_rank, comm_size, proc_name, len(local_data), fps_path))

    if comm_rank == root:
        with utils.timing("rank {}: mol to fps:".format(comm_rank)):
            data_pipeline.mol2fps_mpi(save_path=local_data, fps_path=fps_path)
            sys.stdout.write("process {} done\n".format(comm_rank))
    else:
        data_pipeline.mol2fps_mpi(save_path=local_data, fps_path=fps_path)


if __name__ == "__main__":
    main()
