import tarfile
import os
from loguru import logger
from tqdm import tqdm
from openbabel import openbabel
import concurrent.futures
from data.tools import utils
from data.tools import extractor

WORK_DIR = os.path.curdir
OUT_DIR = os.path.join(WORK_DIR, 'output')


def file_convert(inputfile, informat='pdbqt', outformat='smiles'):
    # TODO: Might optimize by passing param id and dir
    id = os.path.splitext(os.path.basename(inputfile))[0]
    dir = os.path.dirname(inputfile)
    conv= openbabel.OBConversion(inputfile,
                                 os.path.join(dir, id + '.' + outformat))
    conv.OpenInAndOutFiles(informat, outformat)
    conv.Convert()
    conv.CloseOutFile()



def main():
    import argparse

    description = 'A simple command-line interface for virtual flow.'
    parser = argparse.ArgumentParser(description=description)
    # TODO: verbose
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Verbose output')
    parser.add_argument('-d','--data_dir', nargs='+', metavar='<path>',
                        required=True, help='Specify input directory',
                        default='.')
    parser.add_argument('-n', '--n_cpu', type=int, default=1,
                        help='assign the number of processes.')
    # group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-m', '--mode', action='store_const', default='DEBUG',
                       help='Mode selection'
                            'DEBUG: for debug.'
                            'RUN: extract, convert, and calculate fingerprint.'
                            'EXTRACT: only extract tar files.')


    args = parser.parse_args()

    """config"""
    data_dir = args.data_dir
    n_cpu = args.n_cpu

    """logging"""
    log_dir = os.path.join(os.path.curdir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, 'debug_{time}.log')
    logger.add(log_file, retention="12:00")

    """mode selection"""
    if args.mode == 'DEBUG':
        DEBUG_DIR = os.path.join(WORK_DIR, 'debug')
        data_dir = os.path.join(DEBUG_DIR, 'example')



    elif args.mode == 'RUN':

        pass

    elif args.mode == 'EXTRACT':

        pass


    """Extract from .tar files"""
    if args.verbose:
        logger.debug('{} - {}'.format(data_dir,len(os.listdir(data_dir))))
    with utils.tmpdir_manager() as tmp_dir:
        extract(data_dir, tmp_dir, n_cpu)



    """Convert small molecule file format"""
    for d in os.listdir(tmp_dir):




if __name__ == "__main__":
    main()
