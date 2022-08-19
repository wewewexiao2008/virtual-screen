import os
from loguru import logger
from data.tools import utils, pipeline
import mpi4py.MPI as MPI
import sys
import glob

def split_n(origin_list, n):
    res = [[] for i in range(n)]
    for i,k in enumerate(origin_list):
        res[i%n].append(k)
    for i in res:
        yield i

def main():
    import argparse

    description = 'fingerprinting from pdbqt files.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--tmp_dir', type=str, required=True, help='input directory')
    parser.add_argument('-o', '--out_dir', type=str, help='fps output dir', default='./out/')

    args = parser.parse_args()

    out_dir = args.out_dir
    tmp_dir = args.tmp_dir  # './tmp/tmp6sj_2ixi'

    data_pipeline = pipeline.DataPipeline()

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    if comm_rank == 0:
        log_dir = os.path.join(os.path.curdir, 'log')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'debug_{time}.log')
        logger.add(log_file)

        # manager process
        paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)
        send_buf = [i for i in split_n(paths, comm_size)]
        sys.stdout.write("mol num:{}\n".format(len(paths)))
        logger.info("mol num:{}".format(len(paths)))
    else:
        send_buf = None

    local_data = comm.scatter(send_buf, root=0)

    fps_path = os.path.join(out_dir, '{}_{} of {}.fps'.format(proc_name, comm_rank, comm_size))

    with open(fps_path, 'w') as wf:
        wf.writelines("id\tbase64\n")

    sys.stdout.write("process {} of {} on {}, handling {} mols, to {}\n".format(
        comm_rank, comm_size, proc_name, len(local_data), fps_path))

    data_pipeline.mol2fps_mpi(mol_paths=local_data, fps_path=fps_path, rank = comm_rank)

    sys.stdout.write("process {} done".format(comm_rank))


if __name__ == "__main__":
    main()
