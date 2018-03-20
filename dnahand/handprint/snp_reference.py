#!/usr/bin/env python3
from utils import run
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