#!/usr/bin/env python3
from utils import run
import pickle

# def create_pickle(vcf_path, bcftools_bin):
#     snp_reference = {}
#     out = run((f'{bcftools_bin} '
#         f'query -f "%ID\t%REF\t%ALT\t%CHROM\t%POS\n -H " {vcf_path}'))
#     for line in out.split('\n'):
#         if not line:
#             continue
#         ID, REF, ALT, CHR, POS = line.split('\t')
#         snp_reference[ID] = {
#             'REF': REF,
#             'ALT': ALT,
#             'CHR': CHR,
#             'POS': POS,
#         }
#     return snp_reference


def create_snp_data_textfile(reference_vcf_path, snp_data_path, bcftools_bin):
    run((f'{bcftools_bin} query '
        f'-f "%ID\t%REF\t%ALT\t%CHROM\t%POS\n"'
        f'{reference_vcf_path} -o {snp_data_path}'))


def read_snp_data_from_textfile(snp_data_path):
    snp_data = {}
    with open(snp_data_path) as f:
        for line in f:
            if (not line) or (line[0] == '#'):
                continue
            ID, REF, ALT, CHR, POS = line.split('\t')
            snp_data[ID] = {
                'REF': REF,
                'ALT': ALT,
                'CHR': CHR,
                'POS': POS,
            }
        return snp_data


def save_snp_data_as_pickle(snp_data, snp_data_pickle_dest):
    with open (snp_data_pickle_dest, 'wb') as f:
        pickle.dump(snp_data, f)


def load_pickle(pickle_path):
    with open (pickle_path, 'rb') as f:
        return pickle.load(f)