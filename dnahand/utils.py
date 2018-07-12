import os
import subprocess
import sys
import functools
import pandas as pd

import pdb

FINGERPRINT_METHODS = ['sequenom', 'fluidigm']

FIRST_CHROMOSOME = 1
LAST_CHROMOSOME = 22
BCFTOOLS_BIN = '' # Set by pipeline.py
PLINK_BIN = '' # Set by pipeline.py

# def get_reference_vcfs(reference_vcf):
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
        raise Exception
    return result.stdout.decode()


def run_vcf_from_plex(
        vcf_from_plex_bin, filelist, chromosomes, fingerprint_method,
        snpset, vcf_out
    ):
    command = (f'{vcf_from_plex_bin} --input {filelist} '
        f'--chromosomes {chromosomes} --plex_type {fingerprint_method} '
        f'--snpset {snpset} --vcf {vcf_out}')
    run(command)


def gzip_vcfs(vcfs):
    gzipped = []
    print(vcfs)
    for vcf in vcfs:
        run(f'bgzip {vcf}') 
        gzipped.append(f'{vcf}.gz')
    return gzipped


def index_vcf(vcf):
    run(f'{BCFTOOLS_BIN} index -f {vcf}')


def index_vcfs(vcfs):
    for vcf in vcfs:
        index_vcf(vcf)


def convert_sampleids_to_lowercase_vcf(vcf):
    sampleids = run(f'{BCFTOOLS_BIN} query -l {vcf}').split()
    sampleids = {sampleid:sampleid.lower() for sampleid in sampleids if sampleid}
    sample_file = f'{vcf}.sample_names.txt'
    outvcf = f'{vcf}.2'
    with open(sample_file, 'w') as f:
        for old_id, new_id in sampleids.items():
            assert ' ' not in new_id
            assert ' ' not in old_id
            print(old_id, new_id, file=f)
    command = f'{BCFTOOLS_BIN} reheader --samples {sample_file} {vcf} -o {outvcf}'
    run(command)
    os.remove(sample_file)
    os.rename(outvcf, vcf)


def get_subdirectories(directory):
    return [os.path.join(directory, o) for o in os.listdir(directory)
        if os.path.isdir(os.path.join(directory,o))]


def subset_vcf_by_region(vcfin, vcfout, regions):
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
    run(f'{BCFTOOLS_BIN} view -Oz --regions {regions} {vcfin} -o {vcfout}')


def get_sample_ids_from_vcf(vcf):
    sampleids = run(f'{BCFTOOLS_BIN} query -l {vcf}').split()
    # Get rid of any blank lines.
    return {sampleid for sampleid in sampleids if sampleid}


def get_snp_data_from_vcf(vcf_path):
    snp_data = {}
    attributes = '%ID\\t%REF\\t%ALT\\t%CHROM\\t%POS\\n'
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
        handprint_vcfs, reference_snp_pickle, reference_vcf_path,
        subsetted_reference_vcf
    ):
    snp_data = {}
    handprint_rsids = set()
    for handprint_vcf in handprint_vcfs:
        d = get_snp_data_from_vcf(handprint_vcf)
        snp_data = {**snp_data, **d}
    
    regions = {v['CHR'] + ':' + v['POS'] for v in snp_data.values()}
    assert regions
    print(regions)
    subset_vcf_by_region(reference_vcf_path, subsetted_reference_vcf,
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
#         run(f'{BCFTOOLS_BIN} view -Oz -s ^{intersection} {all_finger_prints_to_date} -o {temp}')
#         shutil.move(temp, all_fingerprints_to_date)
#         run(f'{BCFTOOLS_BIN} index -f {all_finger_prints_to_date}')


@functools.lru_cache(maxsize=None)
def get_number_of_snps_from_vcf(vcf_path):
    return len(get_snp_data_from_vcf(vcf_path))


def sort_by_number_of_snps(vcfs):
    key = lambda x: get_number_of_snps_from_vcf(x)
    return sorted(vcfs, key=key, reverse=True)


def append_step_to_vcf_path(vcf, step):
    i = vcf.rfind('.vcf')
    return vcf[:i+1] + step + vcf[i:]


def filter_duplicate_individuals(subsetted_reference_vcf, handprint_vcfs):
    handprint_vcfs = sort_by_number_of_snps(handprint_vcfs)
    handprint_vcfs_processed = []
    sample_ids_so_far = get_sample_ids_from_vcf(subsetted_reference_vcf)
    for vcf in handprint_vcfs:
        sample_ids = get_sample_ids_from_vcf(vcf)
        overlap = sample_ids_so_far & sample_ids
        if overlap:
            vcf_filtered = append_step_to_vcf_path(vcf, 'filtered')
            # overlap = ','.join(overlap)

            exclude_path = f'{vcf}.exclude.txt'
            with open(exclude_path, 'w') as f:
                f.write('\n'.join(overlap))
            run(f'{BCFTOOLS_BIN} view -O v -S ^{exclude_path} {vcf} '
                f'-o {vcf_filtered}')
            vcf = vcf_filtered
        sample_ids_so_far |= sample_ids
        handprint_vcfs_processed.append(vcf)
    return handprint_vcfs_processed


def merge_vcfs(vcf_paths, vcf_out):
    to_merge = ' '.join(vcf_paths)
    run(f'{BCFTOOLS_BIN} merge {to_merge} -o {vcf_out} -Oz')


def prune_by_rsid(vcf_path, exclude_path):
    i = vcf_path.rfind('.vcf.gz')
    out = vcf_path[:i]
    pruned_vcf_path = f'{out}.pruned.vcf.gz'
    regions = get_regions(vcf_path, exclude_path)
    command = f'{BCFTOOLS_BIN} view {vcf_path} -t ^{regions} -o {pruned_vcf_path} -Oz'
    run(command)

    prepruned_vcf_path = f'{out}.prepruned.vcf.gz'

    os.rename(vcf_path, prepruned_vcf_path)
    os.rename(pruned_vcf_path, vcf_path)


def include_by_region_file(vcf_path, include_path):
    i = vcf_path.rfind('.vcf.gz')
    out = vcf_path[:i]
    pruned_vcf_path = f'{out}.pruned.vcf.gz'
    command = f'{BCFTOOLS_BIN} view {vcf_path} -T {include_path} -o {pruned_vcf_path} -Oz'
    run(command)

    prepruned_vcf_path = f'{out}.prepruned.vcf.gz'

    os.rename(vcf_path, prepruned_vcf_path)
    os.rename(pruned_vcf_path, vcf_path)



def get_LD_prune_list(vcf_path):
    i = vcf_path.rfind('.vcf')
    out = vcf_path[:i]
    command = f'{PLINK_BIN} --indep-pairwise 50 5 0.2 --vcf {vcf_path} --out {out}'
    run(command)
    return f'{out}.prune.out'


def get_MAF_exclude_list(frq_path, out_path, maf_min=0.05):
    df = pd.read_csv(frq_path, sep=r'\s+', engine='python')
    subsetted = df[df['MAF'] <= maf_min] 
    subsetted.to_csv(out_path, index=False, columns=['SNP'], header=False)


def get_regions(vcf_path, rsids_exclude_path):
    with open(rsids_exclude_path) as f:
        rsids = [x.strip() for x in f.readlines()]
    lines = run((f'{BCFTOOLS_BIN} query -f "%ID\\t%CHROM\\t%POS\\n" '
        f'{vcf_path}'))
    regions = []
    for line in lines.split('\n'):
        if not line:
            continue
        ID, CHROM, POS = line.split('\t')
        if ID in rsids:
            regions.append(f'{CHROM}:{POS}')

    return ','.join(regions)
