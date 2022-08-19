import multiprocessing
import os
from loguru import logger
from data.tools import utils, pipeline
from model.cluster import Reducer


def main():
    import argparse

    description = 'A simple command-line interface for extracting.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-d', '--data_dir', type=str, metavar='<data_path>', required=True,
                        help='input directory')
    parser.add_argument("-c", "--ckpt_dir", type=str, default='./checkpoint')
    parser.add_argument('-o', '--out_path', type=str, metavar='<out_path>', required=True,
                        help='output directory')
    parser.add_argument('-n', '--n_cpu', type=int,
                        help='assign the number of processes.')
    parser.add_argument('-t', '--type', type=str, default='morgan',
                        help='fingerprint type: morgan or maccs.')
    parser.add_argument('-x', '--delete', type=int, help='delete temp?')
    # group = parser.add_mutually_exclusive_group(required=True)

    args = parser.parse_args()

    """config"""
    data_dir = args.data_dir
    ckpt_dir = args.ckpt_dir
    out_path = args.out_path
    delete = True if args.delete == 1 else False

    n_cpu = args.n_cpu
    fp_type = args.type

    n_groups_l1 = 10000
    n_groups_l2 = 10000

    # logging
    log_dir = os.path.join(os.path.curdir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file)

    if not os.path.exists('./tmp/'):
        os.makedirs('./tmp/')

    if not os.path.exists('./out/'):
        os.makedirs('./out/')

    """main procedure"""
    with utils.tmpdir_manager('./tmp', delete=delete) as tmp_dir:

        data_pipeline = pipeline.DataPipeline(
            data_dir=data_dir,
            n_cpu=n_cpu,
            fp_type=fp_type
        )

        fps_path = os.path.join('./out/', out_path)

        reducer = Reducer(fps_path=fps_path,
                          checkpoint_path=ckpt_dir,
                          n_clusters=n_groups_l1,
                          batch_size=1000,
                          max_iter=500,
                          init_size=n_groups_l1)

        # logger.info("cpu count: {}/{}".format(n_cpu, multiprocessing.cpu_count()))


        with utils.timing("extracting"):
            data_pipeline.extract(tmp_dir)
        with utils.timing("calculating fingerprint"):
            # pack fingerprint as a data file
            data_pipeline.fingerprint(tmp_dir, fps_path)
        #
        # if not os.path.exists(ckpt_dir):
        #     os.makedirs(ckpt_dir)
        #
        # with utils.timing("layer 1: {}".format(n_groups_l1)):
        #     reducer.run()


if __name__ == "__main__":
    main()
