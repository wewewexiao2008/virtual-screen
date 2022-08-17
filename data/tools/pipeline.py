import os, glob
import tarfile
from loguru import logger
from tqdm import tqdm
from concurrent import futures
from concurrent.futures import wait, ALL_COMPLETED

from rdkit.Chem import AllChem
from rdkit.DataStructs import ExplicitBitVect

from scripts.rdkit2pdbqt import MolFromPDBQTBlock


def _extract_single(d_path, tmp_dir, filename):
    """extract a single 'AABBCC' directory"""
    with tarfile.open(os.path.join(d_path, filename)) as tf:
        names = tf.getnames()
        tf.extractall(tmp_dir)
    for name in names:
        dirname = os.path.dirname(name)
        with tarfile.open(os.path.join(tmp_dir, name)) as tf:
            tf.extractall(os.path.join(tmp_dir, dirname))
        os.remove(os.path.join(tmp_dir, name))


def _get_id_from_path(path):
    basename = os.path.basename(path)
    mol_id = os.path.splitext(basename)[0]
    return mol_id


class DataPipeline:
    """Tar extractor, format converter, fingerprint wrapper"""
    def __init__(self,
                 *,
                 data_dir: str,
                 out_dir: str,
                 fp_type: str = 'morgan',
                 n_cpu: int = 16):
        """
        :param data_dir: raw data directory
        :param out_dir: fingerprint output directory
        :param fp_type: 'morgan' or 'maccs'
        :param n_cpu: number of processes
        """
        self.data_dir = data_dir
        self.out_dir = out_dir
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
                logger.warning("can't read mol from {}".format(mol_id))
                return
            base64 = self._mol_to_fingerprint_base64(mol)
            ext = '.fp'
            with open(os.path.join(self.out_dir, mol_id + ext), 'w') as fp_f:
                fp_f.write(base64)

    def _pack_fingerprint(self, filename, max_record):
        """
        Compact fingerprint to one single file for storage convenience.
        :return: None
        """
        fp_ls = glob.glob(r'{}\*.fp'.format(self.out_dir))
        # step = max_record
        # groups = [fp_ls[i:i + step] for i in range(0, len(fp_ls), step)]
        # for i, g in enumerate(fp_ls):
        logger.info("Packing fingerprints {}".format(len(fp_ls)))
        with open(filename, 'w') as wf:
            wf.writelines("id\tbase64\n")
            for k in tqdm(fp_ls):
                with open(k, 'r') as f:
                    basename = os.path.basename(k)
                    mol_id = os.path.splitext(basename)[0]
                    wf.writelines("{}\t{}\n".format(mol_id, f.read()))

    def extract(self, tmp_dir):
        with futures.ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
            d_ls = os.listdir(self.data_dir)
            logger.info("directory num: {}".format(len(d_ls)))
            tasks = []
            for d in d_ls:
                dpath = os.path.join(self.data_dir, d)
                if self.n_cpu > 1:
                    tasks.extend([executor.submit(_extract_single, dpath, tmp_dir, f) for f in (os.listdir(dpath))])
                else:
                    for f in tqdm(os.listdir(dpath)):
                        _extract_single(dpath, tmp_dir, f)
            if self.n_cpu > 1:
                wait(tasks, return_when=ALL_COMPLETED)
            logger.info("extracting done.")

    def fingerprint(self, tmp_dir, pack_filename='packed.fps'):
        """
        .pdbqt file to base64 fingerprint and pack as .fps
        :param max_record:
        :param pack_filename:
        :param tmp_dir: temp directory containing .pdbqt files
        :return: None
        """
        with futures.ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
            # 一级目录 AAAAAA
            paths = glob.glob(r'{}\**\*.pdbqt'.format(tmp_dir), recursive=True)
            logger.info('mol num: {}'.format(len(paths)))

            if self.n_cpu > 1:
                tasks = [executor.submit(self._pdbqt2fingerprint, path) for path in paths]
                        # if 'Z1692919946_1_T1.pdbqt' not in path]
                wait(tasks, return_when=ALL_COMPLETED)
            else:
                for f in tqdm(paths):
                    self._pdbqt2fingerprint(full_path=f)

            # pack
            if pack_filename is not None:
                self._pack_fingerprint(pack_filename)
                logger.info("fingerprinting done.")

