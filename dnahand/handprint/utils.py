from collections import defaultdict
from csv import DictWriter

def get_signature_2_sangerid(sangerid_2_best):
    signature_2_sangerid = defaultdict(set)
    for sangerid, sangerid_2_best in sangerid_2_best.items():
        signature = tuple(sorted([x.get_rsid() for x in sangerid_2_best]))
        signature_2_sangerid[signature].add(sangerid)
    return signature_2_sangerid


def write_snpset(rsids, snp_ref, path):
    with open(path, 'w') as f:
        fieldnames = ['#SNP_NAME', 'REF_ALLELE', 'ALT_ALLELE', 'CHR', 'POS', 'STRAND']
        writer = DictWriter(f, fieldnames, delimiter='\t')
        writer.writeheader()
        for rsid in rsids:
            data = snp_ref[rsid]
            CHR = data['CHR']
            row = {
                '#SNP_NAME': rsid,
                'REF_ALLELE': data['REF'],
                'ALT_ALLELE': data['ALT'],
                'CHR': f'Chr{CHR}',
                'POS': data['POS'],
                'STRAND': '+',
            }
            for v in row.values():
                assert v
            writer.writerow(row)