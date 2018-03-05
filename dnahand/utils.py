import os
import subprocess
import sys

FINGERPRINT_METHODS = ['sequenom', 'fluidigm']
FIRST_CHROMOSOME = 1
LAST_CHROMOSOME = 22

def get_reference_vcfs(reference_vcf_pattern):
    vcf_paths = []
    for chromosome in range(FIRST_CHROMOSOME, LAST_CHROMOSOME + 1):
        vcf_path = reference_vcf_pattern.format(chromosome)
        if os.path.isfile(vcf_path):
            vcf_paths.append(vcf_path)
    return vcf_paths


def run(command):
    print(f'\nrunning: "{command}"')
    # proc = subprocess.Popen(command, shell=True, stdout=stdout,
    #     stderr=stderr)
    # returncode = proc.wait()
    # if returncode != 0:
    #     print(f'Failed {command} with return code {returncode}. Exiting')
    #     sys.exit(-1)
    # (out, err) = proc.communicate()
    # if err:
    #     print(err)
    # return out

    result = subprocess.run(command, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print(result.stderr.decode('utf-8'))
        sys.exit(result.returncode)

    print(result.stderr.decode('utf-8'))


def run_vcf_from_plex(
        vcf_from_plex_bin, filelist, chromosomes, fingerprint_method,
        snpset, vcf_out, stdout, stedrr
    ):
    command = (f'{vcf_from_plex_bin} --input {filelist}'
        f'--chromosomes {chromosomes} --plex_type {fingerprint_method}'
        f'--snpset {snpset} --vcf {vcf_out}')
    run(command, stdout, stedrr)


def gzip_vcfs(vcfs, bgzip_bin):
    gzipped = []
    for vcf in vcfs:
        run(f'bgzip {vcf}') 
        gzipped.append(f'{vcf}.gz')
    return gzipped


def index_vcfs(vcfs, bcftools_bin):
    for vcf in vcfs:
        run(f'{bcftools} index -f {vcf}')


def get_subdirectories(directory):
    return [os.path.join(directory, o) for o in os.listdir(directory)
        if os.path.isdir(os.path.join(directory,o))]


def get_sampleids_from_vcf(vcf):
    return set(run(f'{bcftools} query -l {vcf}').decode().split())


def get_rsids_from_vcf(vcf):
    out = run(f'{bcftools} query -f "%ID\t%REF\t%ALT\t%POS\n" {vcf_path}')
    return set(run(f'{bcftools} query -l {vcf}').decode().split())




def subset_reference_vcfs_on_handprint_snps(handprint_vcfs, reference_snp_pickle, reference_vcf_pattern):
    
    vcf_path = reference_vcf_pattern.format(chromosome)






def remove_duplicate_sample_ids(handprint_vcfs, reference_vcfs):
    handprint_ids = set()
    for handprint_vcf in handprint_vcfs:
        incoming_ids.update(get_sampleids_from_vcf(handprint_vcf))

    all_finger_prints_to_date_ids = get_sampleids_in_vcfgz(all_finger_prints_to_date)
    intersection = incoming_ids & all_finger_prints_to_date_ids
    if intersection:
        print(f'Warning: incoming fingerprints have an intersection with the sample ids from all to date ({all_finger_prints_to_date}):')
        print('\n'.join(intersection))
        print(f'Keeping the newest to replace the fingerprints from {all_finger_prints_to_date}')
        intersection = ','.join(intersection)
        directory, filename = os.path.split(all_finger_prints_to_date)
        temp = os.path.join(directory, 'temp')
        run(f'{bcftools} view -Oz -s ^{intersection} {all_finger_prints_to_date} -o {temp}')
        shutil.move(temp, all_fingerprints_to_date)
        run(f'{bcftools} index -f {all_finger_prints_to_date}')