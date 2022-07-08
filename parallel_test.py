import os
from concurrent import futures
import time
import contextlib
from loguru import logger
from data.tools import pipeline, utils

@contextlib.contextmanager
def timing(msg: str):
    logger.info('Started {}'.format(msg))
    tic = time.time()
    yield
    toc = time.time()
    logger.info('Finished %s in %.3f seconds'%(msg, toc - tic))


LOG_DIR = os.path.join('.', 'log')
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR)

# LOG_FILE = os.path.join(LOG_DIR, 'debug_.log')

# logger.add(LOG_FILE, retention=2)

"""Test for concurrent.futures"""
if __name__ == '__main__':

    out_dir = './debug/output'
    data_pipeline = pipeline.DataPipeline(
        data_dir='./debug/example',
        n_cpu=16,
        mid_format='pdbqt',
        fp_type='maccs',
        out_dir=out_dir
    )
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    logger.add('./log/debug_{time}.log',retention=2)
        # with futures.ProcessPoolExecutor(1) as executor:
        #     for i in range(5, 10000):
        # tmp_dir = ''
    with utils.tmpdir_manager('.', delete=True) as tmp_dir:
        with timing("timing for extracting"):
            data_pipeline.extract(tmp_dir)
        with timing("timing for fingerprinting"):
            data_pipeline.fingerprint(tmp_dir)
        # future = executor.submit(data_pipeline.extract, tmp_dir)
                # print(future.result())

            # input()