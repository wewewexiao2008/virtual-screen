import os
from loguru import logger
from data.tools import utils, pipeline
import mpi4py.MPI as MPI
import sys
import glob
import shutil
import multiprocessing
import joblib


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

    if comm_rank == root:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        log_dir = os.path.join(os.path.curdir, 'log')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'debug_{time}.log')
        logger.add(log_file)

        # manager process
        with utils.timing("counting pdbqt files"):
            sys.stdout.write("counting pdbqt files...\n")
            paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)

        blks = [blk for blk in split_n(paths, n_blk)]
        with utils.timing("saving list"):
            for blk_id, blk in enumerate(blks):
                with open("paths_blk{}.txt".format(blk_id), 'wb') as wf:
                    joblib.dump(blk, wf, protocol=4)

        data = [split_n(ls, comm_size) for ls in split_n(paths, n_blk)]
        sys.stdout.write("mol num:{}, blk num:{}\n".format(len(paths), n_blk))
        logger.info("mol num:{}".format(len(paths)))
    else:
        data = [None for _ in range(n_blk)]

    for block_id, send_block in enumerate(data):
        local_data = comm.scatter([_ for _ in send_block]
                                  if send_block is not None else send_block, root=root)

        fps_path = os.path.join(out_dir,
                                '{}-{}_block{}.fps'.format(proc_name, comm_rank, block_id))
        with open(fps_path, 'w') as wf:
            wf.writelines("id\tbase64\n")

        sys.stdout.write("process {} of {} on {}, handling {} mols, to {}\n".format(
            comm_rank, comm_size, proc_name, len(local_data), fps_path))

        if comm_rank == root:
            with utils.timing("rank {}: mol to fps, block {}:".format(root,block_id)):
                data_pipeline.mol2fps_mpi(
                    mol_paths=local_data, fps_path=fps_path)
                sys.stdout.write("block: {}, process {} done\n".format(block_id, comm_rank))
                # logger.info("block: {}, process {} done\n".format(block_id, comm_rank))
        else:
            data_pipeline.mol2fps_mpi(
                mol_paths=local_data, fps_path=fps_path)
            # sys.stdout.write("block: {}, process: {} done\n".format(block_id, comm_rank))


if __name__ == "__main__":
    main()
