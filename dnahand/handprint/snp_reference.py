#!/usr/bin/env python3
from handprint.utils import run
import pickle

def create_pickle(vcf_path, bcftools):
    snp_reference = {}
    out = run(f'{bcftools} query -f "%ID\t%REF\t%ALT\t%CHROM\t%POS\n" {vcf_path}')
    for line in out.split('\n'):
        if not line:
            continue
        ID, REF, ALT, CHR, POS = line.split('\t')
        snp_reference[ID] = {
            'REF': REF,
            'ALT': ALT,
            'CHR': CHR,
            'POS': POS,
        }
    return snp_reference


def save_pickle(snp_reference, pickle_dest):
    with open (pickle_dest, 'wb') as f:
        pickle.dump(snp_reference, f)


def load_pickle(pickle_path):
    with open (pickle_path, 'rb') as f:
        return pickle.load(f)


def main():
    snp_reference = create_pickle(
        vcf_path='/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/new_imputation_dan/all_studies.sampleids_cleaned_to_lowercase.bcf.gz',
        bcftools='/software/hgi/pkglocal/bcftools-1.6-htslib-1.6-htslib-plugins-6f2229e0-irods-git-4.2.2-plugin_kerberos-2.0.0/bin-wrap/bcftools'
    )
    save_pickle(snp_reference,
        pickle_dest='/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/new_imputation_dan/snp_data_all_studies.pickle'
    )


if __name__ == '__main__':
    main()