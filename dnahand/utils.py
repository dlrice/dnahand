import os
import subprocess
import sys
import functools

FINGERPRINT_METHODS = ['sequenom', 'fluidigm']
FIRST_CHROMOSOME = 1
LAST_CHROMOSOME = 22
BCFTOOLS_BIN = '' # Set by pipeline.py

# def get_reference_vcfs(reference_vcf_pattern):
#     vcf_paths = []
#     for chromosome in range(FIRST_CHROMOSOME, LAST_CHROMOSOME + 1):
#         vcf_path = reference_vcf_pattern.format(chromosome)
#         if os.path.isfile(vcf_path):
#             vcf_paths.append(vcf_path)
#     return vcf_paths


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
        print(result.stderr.decode(), file=sys.stderr)
        sys.exit(result.returncode)
    return result.stdout.decode()


def run_vcf_from_plex(
        vcf_from_plex_bin, filelist, chromosomes, fingerprint_method,
        snpset, vcf_out
    ):
    command = (f'{vcf_from_plex_bin} --input {filelist}'
        f'--chromosomes {chromosomes} --plex_type {fingerprint_method}'
        f'--snpset {snpset} --vcf {vcf_out}')
    run(command)


def gzip_vcfs(vcfs, bgzip_bin):
    gzipped = []
    for vcf in vcfs:
        run(f'bgzip {vcf}') 
        gzipped.append(f'{vcf}.gz')
    return gzipped


def index_vcfs(vcfs):
    for vcf in vcfs:
        run(f'{BCFTOOLS_BIN} index -f {vcf}')


def get_subdirectories(directory):
    return [os.path.join(directory, o) for o in os.listdir(directory)
        if os.path.isdir(os.path.join(directory,o))]


def subset_vcf_by_region(vcfin, vcfout, regions, subset_vcf_by_region):
    """
    Args:
        vcfin (str): Path to vcf to subset.
        vcfout (str): Path to subsetted vcf.
        regions (list(str)): A list of chrom:position to include.

    Returns:
        None

    Side-effects:
        Subsetted vcfout is saved.
    """
    regions = ','.join(regions)
    run(f'{bcftools} view -Oz --regions {regions} {vcfin} -o {vcfout}')


def get_sample_ids_from_vcf(vcf):
    sampleids = run(f'{BCFTOOLS_BIN} query -l {vcf}').split()
    # Get rid of any blank lines.
    return {sampleid in sampleids if sampleid}


def get_snp_data_from_vcf(vcf_path):
    snp_data = {}
    attributes = '%ID\t%REF\t%ALT\t%CHR\t%POS\n'
    out = run(f'{BCFTOOLS_BIN} query -f "{attributes}" {vcf_path}')
    for line in out.split('\n'):
        if not line:
            continue
        ID, REF, ALT, CHR, POS = line.split('\t')
        snp_data[ID] = {
            'REF': REF,
            'ALT': ALT,
            'CHR': CHR,
            'POS': POS,
        }
    return snp_data


def subset_reference_vcfs_on_handprint_snps(
        handprint_vcfs, reference_snp_pickle, reference_vcf_pattern,
        reference_vcf_subset_directory
    ):
    snp_data = {}
    handprint_rsids = set()
    for handprint_vcf in handprint_vcfs:
        d = get_rsids_from_vcf(handprint_vcf)
        snp_data = {**snp_data, **d}
    chromosomes = {v['CHR'] for v in snp_data.values()}
    for chromosome in chromosomes:
        regions = {v['CHR'] + ':' + v['POS'] for v in snp_data.values()
            if v['CHR'] == chromosome}
        if not regions:
            continue
        reference_vcf_path = reference_vcf_pattern.format(chromosome)
        vcf_subset_path = os.path.join(reference_vcf_subset_directory,
            f'reference{chromosome}.vcf.gz')
        subset_vcf_by_region(reference_vcf_path, vcf_subset_path,
            regions)


def concat_vcfs(vcf_paths, vcf_out):
    vcf_paths = ' '.join(vcf_paths)
    command = f'{BCFTOOLS_BIN} concat {vcf_paths} -Oz -o {vcf_out}'
    run(command)

# def remove_duplicate_sample_ids(vcfs):
#     handprint_sample_ids = set()
#     for handprint_vcf in handprint_vcfs:
#         handprint_sample_ids.update(get_sampleids_from_vcf(handprint_vcf))

#     all_finger_prints_to_date_ids = get_sampleids_in_vcfgz(all_finger_prints_to_date)
#     intersection = incoming_ids & all_finger_prints_to_date_ids
#     if intersection:
#         print(f'Warning: incoming fingerprints have an intersection with the sample ids from all to date ({all_finger_prints_to_date}):')
#         print('\n'.join(intersection))
#         print(f'Keeping the newest to replace the fingerprints from {all_finger_prints_to_date}')
#         intersection = ','.join(intersection)
#         directory, filename = os.path.split(all_finger_prints_to_date)
#         temp = os.path.join(directory, 'temp')
#         run(f'{bcftools} view -Oz -s ^{intersection} {all_finger_prints_to_date} -o {temp}')
#         shutil.move(temp, all_fingerprints_to_date)
#         run(f'{bcftools} index -f {all_finger_prints_to_date}')


@functools.lru_cache(maxsize=None)
def get_number_of_snps_from_vcf(vcf_path):
    return len(get_snp_data_from_vcf(vcf_path))


def sort_by_number_of_snps(vcfs, subset_vcf_by_region):
    key = lambda x: get_number_of_snps_from_vcf(x)
    return sorted(vcfs, key=key, reverse=True)


def append_step_to_vcf_path(vcf, step):
    i = vcf.rfind('.vcf.gz')
    return vcf[:i+1] + step + vcf[i:]


def filter_duplicate_individuals(vcfs):
    vcfs_processed = []
    sample_ids_so_far = get_sample_ids_from_vcf(vcfs[0])
    vcfs_processed.append(vcfs[0])
    for vcf in vcfs[1:]:
        sample_ids = get_sample_ids_from_vcf(vcf)
        overlap = sample_ids_so_far & sample_ids
        if overlap:
            vcf_filtered = append_step_to_vcf_path(vcf, 'filtered')
            overlap = ','.join(overlap)
            run(f'{bcftools} view -Oz -s ^{overlap} {vcf} '
                f'-o {vcf_filtered}')
            vcf = vcf_filtered
        sample_ids_so_far |= sample_ids
        vcfs_processed.append(vcf)
    return vcfs_processed


def merge_vcfs(vcf_paths, vcf_out):
    to_merge = ' '.join(vcf_paths)
    run(f'{bcftools} merge {to_merge} -o {vcf_out} -Oz')