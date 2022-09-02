import contextlib
import glob
import shutil
import tempfile
import time
from typing import Optional
from loguru import logger
# from absl
import math

from rdkit.Chem import Draw

from scripts.rdkit2pdbqt import MolFromPDBQTBlock


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None, delete: bool = True):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        if delete:
            logger.info("Deleting temp dir...")
            shutil.rmtree(tmpdir, ignore_errors=True)
        logger.info("done.")


@contextlib.contextmanager
def timing(msg: str):
    logger.info('Started {} '.format(msg))
    tic = time.time()
    yield
    toc = time.time()
    logger.info('Finished %s in %.3f seconds' % (msg, toc - tic))


def show_raw(df, mol_dir, layer1, layer2, out_dir='.', save=False):
    mols = []
    ids = []
    for mol_id in df[(df['layer1'] == layer1) & (df['layer2'] == layer2)]:
        _ls = glob.glob(r'{}/**/{}.pdbqt'.format(mol_dir, mol_id), recursive=True)
        if len(_ls) == 0:
            continue
        mol_path = _ls[0]
        with open(mol_path, 'r') as mol_f:
            mol = MolFromPDBQTBlock(mol_f.read())
        mols.append(mol)
        ids.append(mol_id)

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=4,
        subImgSize=(200, 200),
        legends=[i for i in ids]
    )
    if save:
        img.save('{}/{}_{}.png'.format(out_dir, layer1, layer2))

    return img


def split_n(origin_list, n):
    l = len(origin_list)
    res = [origin_list[math.floor(i / n * l):math.floor((i + 1) / n * l)] for i in range(n)]
    for i in res:
        yield i
