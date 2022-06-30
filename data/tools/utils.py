import contextlib
import shutil
import tempfile
import time
from typing import Optional
from loguru import logger
# from absl


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


@contextlib.contextmanager
def timing(msg: str):
    logger.info('Started {}'.format(msg))
    tic = time.time()
    yield
    toc = time.time()
    logger.info('Finished %s in %.3f seconds'%(msg, toc - tic))
