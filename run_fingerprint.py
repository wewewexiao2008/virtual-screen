import os
from loguru import logger
from data.tools import utils, pipeline
import mpi4py.MPI as MPI
import sys
import glob


def main():
    import argparse

    description = 'fingerprinting from pdbqt files.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-i', '--tmp_dir', type=str, required=True, help='input directory')
    parser.add_argument('-o', '--out_dir', type=str, help='fps output dir', default='./out/')

    args = parser.parse_args()

    out_dir = args.out_dir
    tmp_dir = args.tmp_dir  # './tmp/tmp6sj_2ixi'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    log_dir = os.path.join(os.path.curdir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file)

    data_pipeline = pipeline.DataPipeline()

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    if comm_rank == 0:
        # manager process
        paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)
        sys.stdout.write("mol num:{}".format(len(paths)))
        logger.info("mol num:{}".format(len(paths)))
    else:
        paths = None

    local_data = comm.scatter(paths, root=0)

    fps_path = os.path.join(out_dir, '{}_{} of {}.fps'.format(proc_name, comm_rank, comm_size))

    with open(fps_path, 'w') as wf:
        wf.writelines("id\tbase64\n")

    sys.stdout.write("process {} of {} on {}, handling {} mols, to %s\n".format(
        comm_rank, comm_size, proc_name, len(local_data), fps_path))

    data_pipeline.mol2fps_mpi(mol_paths=local_data, fps_path=fps_path)


if __name__ == "__main__":
    main()
