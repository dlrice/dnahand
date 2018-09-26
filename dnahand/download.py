import subprocess
import os
import utils
import hashlib
import json
from getpass import getpass

password = None

class Sanger_Sample_List_iRODS_DB(object):
    def __init__(self, db_path):
        self.db_path = db_path
        if os.path.exists(db_path):
            with open(db_path) as f:
                self.db = json.load(f)
        else:
             self.db = {}

    def is_sanger_sample_list_seen(self, sanger_sample_list_path):
        path_hash = self.get_sha512_for_file(sanger_sample_list_path)
        return path_hash in self.db.values()

    def set_sanger_sample_list_seen(self, sanger_sample_list_path):
        path_hash = self.get_sha512_for_file(sanger_sample_list_path)
        filename = os.path.basename(sanger_sample_list_path)
        self.db[filename] = path_hash

    def get_sha512_for_file(self, filepath):
        sha512 = hashlib.sha512()
        with open(filepath, 'rb') as f:
            sha512.update(f.read())
        return sha512.hexdigest()

    def close_db(self):
        with open(self.db_path, 'w') as f:
            json.dump(self.db, f, indent=4)


def digest_sample_lists_directory(
        directory, sample_list_irods_db_path, fingerprints_directory,
        baton_bin, baton_metaquery_bin,
        baton_get_bin, irods_credentials_path, n_max_processes=25):
    """
    1. Get the directory of files that list the sample ids to download.
    2. For all files that have not been seen before, download and update shelve
    3.
    """
    print(f'Looking in {directory} to see what needs to be downloaded but skipping entries found in {sample_list_irods_db_path}.')
    db = Sanger_Sample_List_iRODS_DB(sample_list_irods_db_path)

    for fingerprint_method in utils.FINGERPRINT_METHODS:
        fingerprint_method_directory = os.path.join(
            fingerprints_directory, fingerprint_method)
        os.makedirs(fingerprint_method_directory, exist_ok=True)

    dir_contents = os.listdir(directory)
    for file in dir_contents:
        file = os.path.join(directory, file)
        if not os.path.isfile(file):
            print(file, ' is not a file')
            continue
        if not db.is_sanger_sample_list_seen(file):
            for fingerprint_method in utils.FINGERPRINT_METHODS:
                fingerprint_method_directory = os.path.join(
                    fingerprints_directory, fingerprint_method)
                download_fingerprints(file,
                    fingerprint_method_directory, fingerprint_method,
                    baton_bin, baton_metaquery_bin, baton_get_bin,
                    irods_credentials_path, n_max_processes)
            db.set_sanger_sample_list_seen(file)
        else:
            print(f'...{file} already sample IDs already downloaded.')

    db.close_db()


def get_all_sangerids_from_sample_lists_directory(directory):
    sangerids = set()
    dir_contents = os.listdir(directory)
    for file in dir_contents:
        file = os.path.join(directory, file)
        if not os.path.isfile(file):
            print(file, ' is not a file')
            continue
        sangerids.update(read_sangerids(file))
    return sangerids


def read_sangerids(sample_list_path):
    with open(sample_list_path) as f:
        return {x.strip() for x in f.readlines() if x}


def download_fingerprints(
        sample_list_path, fingerprints_directory, fingerprint_method,
        baton_bin, baton_metaquery_bin, baton_get_bin,
        irods_credentials_path, n_max_processes=25,
    ):
    """
    Args:
        sample_list_path (str): Path to a headerless text file which
            lists the Sanger sample IDs to download.
        out_directory (str): Directory where the fingerprints will be
            saved. If already exists returns error.
        fingerprint_method (str, optional): Fingerprint type either:
            'sequenome', 'fluidigm', or 'both'. Defaults to 'both'.
        irods_credentials_path (str, optional): Path to a text file
            containing user's irods password. If not supplied, user will
            be prompted.
        n_max_processes (int, optional): Maximum number of processes to
            use while downloading. The bottleneck when downloading is
            simultaneous irods transfers. Defaults to 25.
    Returns:
        bool: True if successful, False otherwise.

        Creates directory strucutre:
                fingerprints_directory/fingerprint_method: downloaded CSV files
    """

    # Get Sanger Sample IDs
    sangerids = read_sangerids(sample_list_path)
    print(f'...downloading  {len(sangerids)} Sanger IDs from {sample_list_path}')

    login_to_irods(irods_credentials_path)

    # fingerprints_directory = os.path.join(fingerprints_directory, fingerprint_method)
    os.makedirs(fingerprints_directory, exist_ok=True)
    os.chdir(fingerprints_directory)
    processes = set()

    print('...downloading fingerprints from irods')
    for sangerid in sangerids:
        command = (f'{baton_bin} -a sample -v {sangerid} -o = -a {fingerprint_method}_plate -v % -o like |'
            f'{baton_metaquery_bin} -z seq | jq \'.[]\' | {baton_get_bin} --save')
        # print(command)
        processes.add(subprocess.Popen(command, shell=True))
        if len(processes) >= n_max_processes:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

    for p in processes:
       p.wait()


def login_to_irods(path=None):
    global password
    if not password:
        if path:
            with open(path) as f:
                password = f.readline().strip()
        else:
            password = getpass('enter kinit password: ')
    os.system('kinit <<< {}'.format(password))