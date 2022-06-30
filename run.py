import tarfile
import os
from loguru import logger
from tqdm import tqdm
from openbabel import openbabel
import concurrent.futures


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

def extract(data_dir, tmp_dir, num_proc):

    for d in os.listdir(data_dir):
        d_path = os.path.join(data_dir, d)
        for f in tqdm(os.listdir(d_path)):
            if num_proc == 1:
                _extract(d_path, tmp_dir, f)
            elif num_proc > 1:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    executor.map(_extract,)




def _extract(d_path, tmp_dir, f):
    with tarfile.open(os.path.join(d_path, f)) as tf:
        names = tf.getnames()
        tf.extractall(tmp_dir)
    for name in names:
        dirname = os.path.dirname(name)
        with tarfile.open(os.path.join(tmp_dir, name)) as tf:
            tf.extractall(os.path.join(tmp_dir, dirname))


def main():
    import argparse

    description = 'A simple command-line interface for virtual flow.'
    parser = argparse.ArgumentParser(description=description)
    # TODO: verbose
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Verbose output')
    parser.add_argument('-i', '--input', nargs='+', metavar='<path>',
                        required=True, help='Specify input directory',
                        default='.')
    parser.add_argument('-n', '--num_proc', type=int, default=1,
                        help='assign the number of processes.')
    # group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-m', '--mode', action='store_const', default='DEBUG',
                       help='Mode selection'
                            'DEBUG: for debug.'
                            'RUN: extract, convert, and calculate fingerprint.'
                            'EXTRACT: only extract tar files.')


    args = parser.parse_args()

    """config"""
    data_dir = args.input
    tmp_dir = os.path.join(WORK_DIR, 'tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    num_proc = args.num_proc

    """mode selection"""
    if args.mode == 'DEBUG':
        DEBUG_DIR = os.path.join(WORK_DIR, 'debug')
        data_dir = os.path.join(DEBUG_DIR, 'example')

        """logging"""
        log_dir = os.path.join(WORK_DIR, 'log')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        log_file = os.path.join(log_dir, 'debug_{time}.log')
        logger.add(log_file, retention="12:00")
        logger.debug('DEBUG mode, data directory: {}'.format(DEBUG_DIR))

    elif args.mode == 'RUN':
        data_dir = args.input

        pass

    elif args.mode == 'EXTRACT':

        pass


    """Extract from .tar files"""
    if args.verbose:
        logger.debug('{} - {}'.format(data_dir,len(os.listdir(data_dir))))
    extract(data_dir, tmp_dir, num_proc)


    """Convert small molecule file format"""
    for d in os.listdir(tmp_dir):




if __name__ == "__main__":
    main()
