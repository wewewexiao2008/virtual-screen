import argparse
import glob

from loguru import logger
import os
import sys
from sklearn.datasets import make_blobs
from joblib import Parallel, parallel_backend
from joblib import register_parallel_backend
from joblib import delayed
from joblib import cpu_count
from ipyparallel import Client
from ipyparallel.joblib import IPythonParallelBackend
import numpy as np
import datetime
# module in the same directory
from model.cluster import Reducer
from data.tools.utils import timing


FILE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(FILE_DIR)

# prepare the logger
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fps_dir', type=str, required=True)
parser.add_argument('-o', '--out_dir', type=str, required=True)
parser.add_argument("-p", "--profile", default="ipy_profile",
                    help="Name of IPython profile to use")

args = parser.parse_args()
profile = args.profile
fps_dir = args.fps_dir
out_dir = args.out_dir

logger.add(os.path.join(FILE_DIR, profile + '.log'))
logger.info("number of CPUs found: {0}".format(cpu_count()))
logger.info("args.profile: {0}".format(profile))

# prepare the engines
c = Client(profile=profile)
NB_WORKERS = int(os.environ.get("NB_WORKERS", 1))
# wait for the engines
c.wait(timeout=5)

# The following command will make sure that each engine is running in
# the right working directory to access the custom function(s).
c[:].map(os.chdir, [FILE_DIR] * len(c))
logger.info("c.ids :{0}".format(str(c.ids)))
bview = c.load_balanced_view()
register_parallel_backend('ipyparallel',
                          lambda: IPythonParallelBackend(view=bview))


# prepare it for the custom function
# some parameters to test in parallel
param_space = {
    'paths': glob.glob(r'{}/*.fps'.format(fps_dir), recursive=True)

}

nc_layer1 = 576 * 2
nc_layer2 = 576 * 4

l1_reducer = Reducer(
    n_clusters=nc_layer1,
    batch_size=1000,
    max_iter=1000,
    init_size=nc_layer1
)
l2_reducer = Reducer(
    n_clusters=nc_layer2,
    batch_size=1000,
    max_iter=1000,
    init_size=nc_layer2
)

# with timing("parallel backend"):
with parallel_backend('ipyparallel'):
    inertia = Parallel(n_jobs=len(c))(
        delayed(l1_reducer.run_with_fps)(fps_path, True)
        for fps_path in param_space['paths'])

logger.info("done")
# write down the number of clusters and the total inertia in a file.
with open(os.path.join(FILE_DIR, 'scores_kmeans.csv'), 'w') as f:
    f.write('fps_path,inertia,\n')
    # f.write("")
    f.write("\n".join(','.join(str(c) for c in _) for _ in inertia))
    f.write('\n')
c.shutdown()
