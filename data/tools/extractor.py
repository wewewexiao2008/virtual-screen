import os
from concurrent import futures

from tqdm import tqdm


class Extractor:
    """Tar files extractor"""

    def __init__(self,
                 *,
                 data_dir,
                 tmp_dir, ):

    def extract(self, data_dir, tmp_dir, num_proc):

        for d in os.listdir(data_dir):
            d_path = os.path.join(data_dir, d)
            for f in tqdm(os.listdir(d_path)):
                if num_proc == 1:
                    self._extract(d_path, tmp_dir, f)
                elif num_proc > 1:
                    with futures.ProcessPoolExecutor() as executor:
                        executor.map(self._extract,)

    def _extract(slef, d_path, tmp_dir, f):
        with tarfile.open(os.path.join(d_path, f)) as tf:
            names = tf.getnames()
            tf.extractall(tmp_dir)
        for name in names:
            dirname = os.path.dirname(name)
            with tarfile.open(os.path.join(tmp_dir, name)) as tf:
                tf.extractall(os.path.join(tmp_dir, dirname))