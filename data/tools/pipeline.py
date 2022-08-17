import os, glob
import tarfile
from loguru import logger
from tqdm import tqdm
from concurrent import futures
from concurrent.futures import wait, ALL_COMPLETED
import multiprocessing

from rdkit.Chem import AllChem
from rdkit.DataStructs import ExplicitBitVect

from scripts.rdkit2pdbqt import MolFromPDBQTBlock


def _extract_single(tar_path, tmp_dir):
    """extract a single 'AABBCC' directory"""
    with tarfile.open(tar_path) as tf:
        names = tf.getnames()
        tf.extractall(tmp_dir)
    for name in names:
        dirname = os.path.dirname(name)
        with tarfile.open(os.path.join(tmp_dir, name)) as tf:
            tf.extractall(os.path.join(tmp_dir, dirname))
        os.remove(os.path.join(tmp_dir, name))
    # logger.info("{}".format(tar_path))


def _get_id_from_path(path):
    basename = os.path.basename(path)
    mol_id = os.path.splitext(basename)[0]
    return mol_id


class DataPipeline:
    """Tar extractor, format converter, fingerprint wrapper"""
    def __init__(self,
                 *,
                 data_dir: str,
                 fp_type: str = 'morgan',
                 n_cpu: int):
        """
        :param data_dir: raw data directory
        :param fp_type: 'morgan' or 'maccs'
        :param n_cpu: number of processes
        """
        self.data_dir = data_dir
        self.fp_type = fp_type
        self.n_cpu = n_cpu
        self.fps = []

    def _mol_to_fingerprint_base64(self, mol):
        if self.fp_type == 'morgan':
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2)
        elif self.fp_type == 'maccs':
            fp = AllChem.GetMACCSKeysFingerprint(mol)
        else:
            fp = ExplicitBitVect(0)
            logger.error("get wrong fingerprint type as {}!".format(self.fp_type))
        base64 = ExplicitBitVect.ToBase64(fp)
        return base64

    def _pdbqt2fingerprint(self, full_path):
        """
        Calculate .pdbqt file's fingerprint.
        Wrapper for parallel processing
        :param full_path:
        :param mid_format: 'sdf' or 'smiles'
        :return: None
        """
        # directly read as rdkit.Chem.rdchem.Mol
        with open(full_path) as mol_f:
            mol_id = _get_id_from_path(full_path)
            mol = MolFromPDBQTBlock(mol_f.read())
            if mol is None:
                # logger.warning("can't read mol from {}".format(mol_id))
                return None
            base64 = self._mol_to_fingerprint_base64(mol)
            return mol_id, base64

    def extract(self, tmp_dir):
        # with futures.ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
        d_ls = os.listdir(self.data_dir)
        paths = glob.glob(r'{}/**/*.tar'.format(self.data_dir), recursive=True)
        logger.info("directory num: {}".format(len(d_ls)))
        logger.info('tarfile num: {}'.format(len(paths)))

        pool = multiprocessing.Pool(self.n_cpu)
        for path in paths:
            pool.apply_async(func=_extract_single, args=(path, tmp_dir, ))
        pool.close()
        pool.join()

    def fingerprint(self, tmp_dir, pack_filename='packed.fps'):
        """
        .pdbqt file to base64 fingerprint and pack as .fps
        :param pack_filename:
        :param tmp_dir: temp directory containing .pdbqt files
        :return: None
        """
        def write_callback(res):
            if res is not None:
                with open(pack_filename, 'a') as wf:
                    line = "{}\t{}\n".format(res[0], res[1])
                    wf.write(line)

        # with futures.ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
            # 一级目录 AAAAAA
        paths = glob.glob(r'{}/**/*.pdbqt'.format(tmp_dir), recursive=True)
        logger.info('mol num: {}'.format(len(paths)))

        with open(pack_filename, 'w') as wf:
            wf.writelines("id\tbase64\n")

        pool = multiprocessing.Pool(self.n_cpu)
        for path in paths:
            pool.apply_async(func=self._pdbqt2fingerprint, args=(path,), callback=write_callback)
        pool.close()
        pool.join()


