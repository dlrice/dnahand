from glob import glob
import os
from collections import defaultdict
import sys
from csv import DictReader, DictWriter
import snp_reference
from utils import *


def get_fieldnames():
    return ['assay', 'snp_assayed', 'x_allele', 'y_allele', 'sample_name',
            'type', 'auto', 'confidence', 'final', 'converted_call', 
            'x_intensity', 'y_intensity']


class FingerprintRow(object):
    def __init__(self, row):
        self.row = row

    def get_assay(self):
        return self.row['assay']
    
    def set_assay(self, assay):
        self.row['assay'] = assay
    
    def reset_assay(self):
        assay = self.get_assay()
        _, a = assay.split('-')
        assay = f'S000-{a}'
        self.set_assay(assay)
    
    def get_row(self):
        return self.row
    
    def get_sangerid(self):
        return self.row['sample_name']
    
    def get_rsid(self):
        return self.row['snp_assayed']
    
    def get_alleles(self):
        return [self.row['x_allele'], self.row['y_allele']]

    def has_alleles(self, alleles):
        return set(self.get_alleles()) == set(alleles)

    def change_strand(self):
        x_allele = self.row['x_allele']
        self.row['x_allele'] = BASE_PAIR_RULE[x_allele]
        y_allele = self.row['y_allele']
        self.row['y_allele'] = BASE_PAIR_RULE[y_allele]

    def get_confidence(self):
        return self.row['confidence']
    
    def __repr__(self):
        sangerid = self.get_sangerid()
        alleles = self.get_alleles()
        confidence = self.get_confidence()
        return f'{sangerid} - {alleles} {confidence}'
    
    
def get_sangerid_2_rsid_2_fingerprintrows(
        fingerprints, fieldnames, snp_ref
    ):
    sangerid_2_rsid_2_fingerprintrows = defaultdict(lambda: defaultdict(list))
    for fingerprint in fingerprints:
        with open(fingerprint) as f:
            reader = DictReader(f, fieldnames, delimiter='\t')
            for row in reader:
                fr = FingerprintRow(row)
                rsid = fr.get_rsid()
                sangerid = fr.get_sangerid()
                if rsid in snp_ref:
                    t = snp_ref[rsid]
                    ref_alleles = [t['REF'], t['ALT']]
                    if fr.has_alleles(ref_alleles):
                        sangerid_2_rsid_2_fingerprintrows[sangerid][rsid].append(fr)
                    else:
                        fr.change_strand()
                        if fr.has_alleles(ref_alleles):
                            sangerid_2_rsid_2_fingerprintrows[sangerid][rsid].append(fr)
                    
    return sangerid_2_rsid_2_fingerprintrows


def get_sangerid_2_bestfingerprintrows(
        sangerid_2_rsid_2_fingerprintrows, snp_ref
    ):
    sangerid_2_bestfingerprintrows = defaultdict(list)
    for sangerid, rsid_2_fingerprintrows in sangerid_2_rsid_2_fingerprintrows.items():
        for rsid, fingerprintrows in rsid_2_fingerprintrows.items():
                frs = sorted(fingerprintrows, key=lambda x: x.get_confidence(), reverse=True)
                sangerid_2_bestfingerprintrows[sangerid].append(frs[0])
    return sangerid_2_bestfingerprintrows


def write_fingerprint(fingerprint_rows, fieldnames, path):
    [x.reset_assay() for x in fingerprint_rows]
    fingerprint_rows = sorted(fingerprint_rows, key=lambda x: x.get_assay())
    with open(path, 'w') as f:
        writer = DictWriter(f, fieldnames, delimiter='\t')
        for fingerprint_row in fingerprint_rows:
            writer.writerow(fingerprint_row.row)


def write_best_merged_fingerprints(
        signature_2_sangerid, sangerid_2_bestfingerprintrows, fieldnames,
        out_directory, snp_ref
    ):
    for index, (signature, sangerids) in enumerate(signature_2_sangerid.items()):
        path = os.path.join(out_directory, f'handprint{index}')
        os.makedirs(path, exist_ok=True)
        snpset_path = os.path.join(path, 'snpset.tsv')
        write_snpset(signature, snp_ref, snpset_path)
        paths = []
        for sangerid in sangerids:
            fingerprint_path = os.path.join(path, f'{sangerid}.csv')
            bestfingerprintrows = sangerid_2_bestfingerprintrows[sangerid]
            write_fingerprint(bestfingerprintrows, fieldnames, fingerprint_path)
            paths.append(str(fingerprint_path))
        filelist = os.path.join(path, f'filelist.txt')
        with open(filelist, 'w') as f:
            f.write('\n'.join(paths))


def generate(
        fingerprints_directory,
        out_directory, snp_reference_pickle
    ):
    """
    Collate and prioritize SNPs across fingerprints
        Args:
            fingerprints_directory (str): Directory which contains all of
                the fingerprints.
            out_directory (str): Directory where the collated fingerprints
                will be saved. If already exists returns error.
            snp_reference_pickle (str): Path to a pickled dictionary of snp
                reference data. The SNP Set is made from SNPs present in 
                both the fingerprints and this pickle. Structure:
                    rsid (str): {
                        'REF': (str),
                        'ALT': (str),
                        'POS': (int),
                        'CHR': (str),
                    }
        Returns:
            bool: True if successful, False otherwise.
    """
    snp_ref = snp_reference.load_pickle(snp_reference_pickle)
    fingerprints = glob(os.path.join(fingerprints_directory, '*.csv'))
    n_fingerprints = len(fingerprints)
    print('Found', n_fingerprints, 'fluidigm fingerprints.')
    if not n_fingerprints:
        print('...SKIPPING')
        return False
    
    fieldnames = get_fieldnames()

    sangerid_2_rsid_2_fingerprintrows = get_sangerid_2_rsid_2_fingerprintrows(
        fingerprints, fieldnames, snp_ref)

    n_sangerids = len(sangerid_2_rsid_2_fingerprintrows)
    print(f'Collected {n_sangerids} sangerids')

    sangerid_2_bestfingerprintrows = get_sangerid_2_bestfingerprintrows(
        sangerid_2_rsid_2_fingerprintrows, snp_ref)

    signature_2_sangerid = get_signature_2_sangerid(
        sangerid_2_bestfingerprintrows)

    print('Found', len(signature_2_sangerid), 'unique snpsets')

    write_best_merged_fingerprints(
        signature_2_sangerid, sangerid_2_bestfingerprintrows, 
        fieldnames, out_directory, snp_ref)

    return True


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Takes directory of Fluidigm \
#         fingerprints, groups by sanger id, selects highest quality rows, \
#         and creates snpsets for each collection of sanger ids which share \
#         the snps after filtering.')
#     parser.add_argument('--fluidigm-dir', help='Directory of fluidigm \
#         fingerprints (csvs).')
#     parser.add_argument('--snp-reference-pickle', help='Path to the pickle \
#         that contains dictionary with rsids: {REF,ALT,POS,chrom} for the \
#         reference collection of snps.')
#     parser.add_argument('--dest-dir', help='Destination for results.')
#     args = vars(parser.parse_args())

#     print('Got arguments:')
#     for k, v in args.items():
#         print(k,v)

#     main(**args)