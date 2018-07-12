import subprocess
import os

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
    print(f'...found {len(sangerids)} Sanger IDs')

    login_to_irods(irods_credentials_path)

    fingerprints_directory = os.path.join(fingerprints_directory, fingerprint_method)
    os.makedirs(fingerprints_directory)
    os.chdir(fingerprints_directory)
    processes = set()

    print('...downloading fingerprints from irods')
    for sangerid in sangerids:
        command = (f'{baton_bin} -a sample -v {sangerid} -o = -a {fingerprint_method}_plate -v % -o like |'
            f'{baton_metaquery_bin} -z seq | jq \'.[]\' | {baton_get_bin} --save')
        print(command)
        processes.add(subprocess.Popen(command, shell=True))
        if len(processes) >= n_max_processes:
            os.wait()
            processes.difference_update([p for p in processes if p.poll() is not None])

    for p in processes:
       p.wait()
    

def login_to_irods(path=None):
    if path:
        with open(path) as f:
            password = f.readline().strip()
        os.system('kinit <<< {}'.format(password))
    else:
        os.system('kinit <<< {}'.format(getpass('enter kinit password: ')))