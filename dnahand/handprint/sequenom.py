from glob import glob
import os
from collections import defaultdict
import sys
from csv import DictReader, DictWriter
import handprint.snp_reference as snp_reference
from handprint.utils import get_signature_2_sangerid, write_snpset
from utils import *

BASE_PAIR_RULE = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}

class FingerprintRowPair(object):
    def __init__(self, fingerprint_rows):
        if len(fingerprint_rows) != 2:
            raise Exception('Only using biallelic SNPs')
        self.fingerprint_rows = fingerprint_rows
        
        s1 = self.fingerprint_rows[0].get_sangerid()
        s2 = self.fingerprint_rows[1].get_sangerid()
        assert s1 == s2
        r1 = self.fingerprint_rows[0].get_rsid()
        r2 = self.fingerprint_rows[1].get_rsid()
        assert r1 == r2

    def get_sangerid(self):
        return self.sangerid
                
    def get_rows(self):
        return [x.get_row() for x in self.fingerprint_rows]
    
    def has_alleles(self, alleles):
        a1 = self.fingerprint_rows[0].get_allele()
        a2 = self.fingerprint_rows[1].get_allele()
        return set([a1, a2]) == alleles
    
    def change_strand(self):
        self.fingerprint_rows[0].change_strand()
        self.fingerprint_rows[1].change_strand()
        
    def get_quality(self):
        q1 = self.fingerprint_rows[0].get_quality()
        q2 = self.fingerprint_rows[1].get_quality()
        return q1 + q2 / 2
        
    def get_rsid(self):
        return self.fingerprint_rows[0].get_rsid()
    
    def __repr__(self):
        sangerid = self.fingerprint_rows[0].get_sangerid()
        a1 = self.fingerprint_rows[0].get_allele()
        a2 = self.fingerprint_rows[1].get_allele()
        q = self.get_quality()
        return f'{sangerid}: {a1} {a2} | {q}'


class FingerprintRow(object):
    def __init__(self, row):
        row['ASSAY_ID'] = row['ASSAY_ID'].split('-')[1]
        self.row = row

    def get_row(self):
        return self.row
    
    def get_sangerid(self):
        return self.row['SAMPLE_ID']
    
    def get_rsid(self):
        return self.row['ASSAY_ID']
    
    def get_plate(self):
        return self.row['PLATE']
    
    def get_allele(self):
        return self.row['ALLELE']
    
    def set_allele(self, allele):
        self.row['ALLELE'] = allele
        
    def change_strand(self):
        allele = self.get_allele()
        self.row['ALLELE'] = BASE_PAIR_RULE[allele]

    def get_status(self):
        return self.row['STATUS']
    
    def get_quality(self):
        """
        The best achievable is 0 (A.Conservative)
        The worst achievable is 25 (Z. ... I am not sure what this is)
        
        A.Conservative
        B.Moderate
        C.Aggressive
        D.Low Probability
        E.User Call
        F.User Call
        I.Bad Spectrum
        N.No-Alleles

        """
        lo = ord('A')
        hi = ord('Z')
        q = ord(self.get_status()[0].upper())
        assert lo <= q <= hi, f'Quality error: {lo} <= {q} <= {hi}'
        return q - lo
        
    def __repr__(self):
        sangerid = self.get_sangerid()
        allele = self.get_allele()
        status = self.get_status()
        return f'{sangerid} - {allele} {status}'

    
def get_sangerid_2_rsid_2_plate_2_fingerprintrows(fingerprints, snp_ref):
    sangerid_2_rsid_2_plate_2_fingerprintrows = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for fingerprint in fingerprints:
        with open(fingerprint) as f:
            reader = DictReader(f, delimiter='\t')
            for row in reader:
                fpr = FingerprintRow(row)
                rsid = fpr.get_rsid()
                sangerid = fpr.get_sangerid()
                plate = fpr.get_plate()
                if rsid in snp_ref:
                       sangerid_2_rsid_2_plate_2_fingerprintrows[sangerid][rsid][plate].append(fpr)

    return sangerid_2_rsid_2_plate_2_fingerprintrows


def get_fieldnames(fingerprint):
    # ASSUMPTION: all fingerprints have the same header
    with open(fingerprint[0]) as f:
        reader = DictReader(f, delimiter='\t')
        return reader.fieldnames


def get_sangerid_2_bestfingerprintrowpairs(sangerid_2_rsid_2_plate_2_fingerprintrows, snp_ref):
    sangerid_2_bestfingerprintrowpairs = defaultdict(list)
    for sangerid, rsid_2_plate_2_fingerprintrows in sangerid_2_rsid_2_plate_2_fingerprintrows.items():
        for rsid, plate_2_fingerprintrows in rsid_2_plate_2_fingerprintrows.items():
            frps = []
            for plate, fingerprintrows in plate_2_fingerprintrows.items():
                if len(fingerprintrows) != 2:
                    continue
                frp = FingerprintRowPair(fingerprintrows)
                t = snp_ref[rsid]
                ref_alleles = set([t['REF'], t['ALT']])
                if frp.has_alleles(ref_alleles):
                    frps.append(frp)
                else:
                    frp.change_strand()
                    if frp.has_alleles(ref_alleles):
                        frps.append(frp)
            if frps:
                frps = sorted(frps, key=lambda x: x.get_quality())
                sangerid_2_bestfingerprintrowpairs[sangerid].append(frps[0])
    return sangerid_2_bestfingerprintrowpairs


def write_fingerprint(fingerprint_row_pairs, fieldnames, path):
    with open(path, 'w') as f:
        writer = DictWriter(f, fieldnames, delimiter='\t')
        writer.writeheader()
        for fingerprint_row_pair in fingerprint_row_pairs:
            for row in fingerprint_row_pair.get_rows():
                ASSAY_ID = row['ASSAY_ID']
                row['ASSAY_ID'] = f'W00000-{ASSAY_ID}'
                writer.writerow(row)


def write_best_merged_fingerprints(signature_2_sangerid, sangerid_2_bestfingerprintrowpairs, fieldnames, out_directory, snp_ref):
    for index, (signature, sangerids) in enumerate(signature_2_sangerid.items()):
        path = os.path.join(out_directory, f'handprint{index}')
        os.makedirs(path, exist_ok=True)
        snpset_path = os.path.join(path, 'snpset.tsv')
        write_snpset(signature, snp_ref, snpset_path)
        paths = []
        for sangerid in sangerids:
            fingerprint_path = os.path.join(path, f'{sangerid}.csv')
            bestfingerprintrowpairs = sangerid_2_bestfingerprintrowpairs[sangerid]
            write_fingerprint(bestfingerprintrowpairs, fieldnames, fingerprint_path)
            paths.append(str(fingerprint_path))
        filelist = os.path.join(path, f'filelist.txt')
        with open(filelist, 'w') as f:
            f.write('\n'.join(paths))


def generate(
        fingerprints_directory, out_directory, snp_reference_pickle
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

# def main(sequenom_dir, snp_reference_pickle, dest_dir):
    
    snp_ref = snp_reference.load_pickle(snp_reference_pickle)
    fingerprints = glob(os.path.join(fingerprints_directory, '*.csv'))
    n_fingerprints = len(fingerprints)
    print('Found', n_fingerprints, 'sequenom fingerprints.')
    if not n_fingerprints:
        print('...SKIPPING')
        return False

    sangerid_2_rsid_2_plate_2_fingerprintrows = get_sangerid_2_rsid_2_plate_2_fingerprintrows(
        fingerprints, snp_ref)

    sangerid_2_bestfingerprintrowpairs = get_sangerid_2_bestfingerprintrowpairs(
        sangerid_2_rsid_2_plate_2_fingerprintrows, snp_ref)

    n_sangerids = len(sangerid_2_bestfingerprintrowpairs)
    print(f'Collected {n_sangerids} sangerids')

    signature_2_sangerid = get_signature_2_sangerid(
        sangerid_2_bestfingerprintrowpairs)

    print(f'Found {len(signature_2_sangerid)} unique snpsets')
    fieldnames = get_fieldnames(fingerprints)
    write_best_merged_fingerprints(signature_2_sangerid, 
        sangerid_2_bestfingerprintrowpairs, fieldnames, out_directory,
        snp_ref)

    return True


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Takes directory of Sequenom \
#         fingerprints, groups by sanger id, selects highest quality rows, \
#         and creates snpsets for each collection of sanger ids which share \
#         the snps after filtering.')
#     parser.add_argument('--sequenom-dir', help='Directory of sequenom \
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