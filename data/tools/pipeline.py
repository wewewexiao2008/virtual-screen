import os,glob
import tarfile
from openbabel import openbabel, pybel
from loguru import logger
from concurrent import futures
from concurrent.futures import  wait, ALL_COMPLETED, FIRST_COMPLETED

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ExplicitBitVect
from tqdm import tqdm

from scripts.rdkit2pdbqt import MolFromPDBQTBlock


def _mol_to_fingerprint_base64(mol, fp_type):
    if fp_type == 'morgan':
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2)
    elif fp_type == 'maccs':
        fp = AllChem.GetMACCSKeysFingerprint(mol)
    else:
        fp = ExplicitBitVect(0)
        logger.error("get wrong fingerprint type as {}!".format(fp_type))
    base64 = ExplicitBitVect.ToBase64(fp)
    return base64


def _convert_single_file(full_path, out_ext):
    """
    Convert file format for a single pdbqt file
    example:
        ./tmp/AABABD/00000/Z18500480_1_T1.pdbqt
    =>  ./tmp/AABABD/00000/Z18500480_1_T1.sdf
    """
    out_full_path = os.path.splitext(full_path)[0] + '.' + out_ext
    # logger.debug(full_path)
    # logger.debug(out_full_path)
    f = pybel.readfile('pdbqt', full_path)
    f.__next__().write(out_ext, out_full_path, overwrite=True)
    f.close()



def _cal_single_fingerprint(path, fp_type, out_dir):
    """
    calculate single fingerprint from file_id
    """
    be = os.path.splitext(path)
    basename = os.path.basename(be[0])
    ext = be[1]
    if ext == '.smiles':
        sup = Chem.SmilesMolSupplier(path)
    elif ext == '.sdf':
        sup = Chem.SDMolSupplier(path)
    else:
        raise TypeError("{}: unknown file format".format(path))

    base64 = [_mol_to_fingerprint_base64(mol, fp_type) for mol in sup]

    with open(os.path.join(out_dir, basename + '.fp'), 'w') as f:
        f.write(base64[0])


def _mol2fingerprint(full_path, ext, out_dir, fp_type="morgan"):
    """
    Wrapper for parallel processing
    :param full_path:
    :param ext: 'sdf' or 'smiles'
    :param out_dir:
    :param fp_type: 'morgan' or 'maccs'
    :return:
    """

    if ext == 'pdbqt':
        with open(full_path) as f:
            mol = MolFromPDBQTBlock(f.read())
            base64 = _mol_to_fingerprint_base64(mol, fp_type)
            basename = os.path.basename(full_path)
            with open(os.path.join(out_dir, basename + '.fp'), 'w') as f:
                f.write(base64[0])
    else:
        _convert_single_file(full_path, ext)
        _path = os.path.splitext(full_path)[0] + '.' + ext
        logger.debug(_path)
        _cal_single_fingerprint(_path, fp_type, out_dir)


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


class DataPipeline:
    """Tar extractor, format converter, fingerprint wrapper"""
    def __init__(self,
                 *,
                 data_dir: str,
                 out_dir: str,
                 n_cpu: int = 1):
        self.data_dir = data_dir
        self.out_dir = out_dir
        self.n_cpu = n_cpu

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

    def fingerprint(self, tmp_dir):
        with futures.ProcessPoolExecutor(max_workers=self.n_cpu) as executor:
            # 一级目录 AAAAAA
            paths = glob.glob(r'{}\**\*.pdbqt'.format(tmp_dir), recursive=True)
            logger.info('mol num: {}'.format(len(paths)))

            if self.n_cpu > 1:
                tasks = [executor.submit(_mol2fingerprint, path, 'pdbqt', self.out_dir) for path in paths]
                        # if 'Z1692919946_1_T1.pdbqt' not in path]
                wait(tasks, return_when=ALL_COMPLETED)
            else:
                for f in tqdm(paths):
                    _mol2fingerprint(full_path=f,
                                     ext='sdf',
                                     out_dir=self.out_dir)
            logger.info("fingerprinting done.")
            out_files = glob.glob(r'{}\*.fp'.format(self.out_dir))
            logger.info('output num: {}'.format(len(out_files)))




