import os
import time
import contextlib
from loguru import logger
from data.tools import utils
from data import pipeline


@contextlib.contextmanager
def timing(msg: str):
    logger.info('Started {}'.format(msg))
    tic = time.time()
    yield
    toc = time.time()
    logger.info('Finished %s in %.3f seconds' % (msg, toc - tic))


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
        fp_type='morgan',
        out_dir=out_dir
    )
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    logger.add('./log/debug_{time}.log', retention=2)
        # with futures.ProcessPoolExecutor(1) as executor:
        #     for i in range(5, 10000):
        # tmp_dir = ''
    with utils.tmpdir_manager('./debug/', delete=True) as tmp_dir:
        with timing("extracting and fingerprinting"):
            data_pipeline.extract(tmp_dir)
            data_pipeline.fingerprint(tmp_dir)
        # future = executor.submit(data_pipeline.extract, tmp_dir)
                # print(future.result())

            # input()