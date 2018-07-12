#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white', color_codes=True)
from collections import defaultdict
from csv import DictReader
from utils import run, get_sample_ids_from_vcf
import os
import pickle


def get_counts_from_csv(kinship_results_csv_path, reference_samples, handprint_samples, dest_dir, duplicate_threshold=0.4):
    KINSHIP_CSV_HEADER = ['pair', 'ibd0', 'ibd1', 'ibd2', 'kinship', 'snps']

    # Holds only the counts of (snps, kinship) pairs.
    nref_2_snp_kinship_2_count = {
        0: defaultdict(int), # All handprint
        1: defaultdict(int), # 1 handprint / 1 reference
        2: defaultdict(int), # All reference
    }

    # Holds the sangerids of those >= duplicate_threshold.
    ref_2_snp_kinship_2_duplicates = {
        0: [], # All handprint
        1: [], # 1 handprint / 1 reference
        2: [], # All reference
    }

    with open(kinship_results_csv_path) as f:
        reader = DictReader(f, delimiter=' ',
            fieldnames=KINSHIP_CSV_HEADER)
        
        for row in reader:
            kinship = row['kinship'][:4]
            snps = row['snps']
            id1, id2, _ = row['pair'].split('\t')
            
            # Sanity check to make sure it belongs to one of the two.
            reference_count = 0
            for _id in (id1, id2):
                if _id in reference_samples:
                    reference_count += 1
                    assert _id not in handprint_samples, \
                        (f'{_id} in reference_samples and '
                        f'handprint_samples')
                else:
                    assert _id in handprint_samples, \
                        (f'{_id} not in reference_samples or '
                        f'handprint_samples')
            
            k = (snps, kinship)
            nref_2_snp_kinship_2_count[reference_count][k] += 1

            float_kinship = float(kinship)
            if float_kinship >= duplicate_threshold:
                _ = (id1, id2, float_kinship, snps)
                ref_2_snp_kinship_2_duplicates[reference_count].append(_)

    path = os.path.join(dest_dir, 'nref_2_snp_kinship_2_count.pickle')
    with open (path, 'wb') as f:
        pickle.dump(nref_2_snp_kinship_2_count, f)


    path = os.path.join(dest_dir, 'ref_2_snp_kinship_2_duplicates.pickle')
    with open (path, 'wb') as f:
        pickle.dump(ref_2_snp_kinship_2_duplicates, f)


if __name__ == '__main__':
    main()

# def get_counts_from_csv(kinship_results_csv_path):
#     counts = defaultdict(int)
#     with open(kinship_results_csv_path) as f:
#         reader = DictReader(f, delimiter=' ',
#             fieldnames=KINSHIP_CSV_HEADER)
#         for row in reader:
#             kinship = row['kinship'][:4]
#             snps = row['snps']
#             counts[(snps, kinship)]+=1

#     return counts


# def numify(x, n):
#     return int(x/n)*n


# def bucket_kinship_nsnps_counts(counts, dnsnp=10, dkinship=0.025):
#     bcounts = defaultdict(int)
#     for k, v in counts.items():
#         snp = int(k[0])
#         snp = numify(snp, dnsnp)
#         kinship = float(k[1])
#         kinship = numify(kinship, dkinship)
#         bcounts[(snp, kinship)] += v
#     return bcounts


# def plot_counts(counts, dest)
#     X, Y, C = [], [], []
#     for (snp, kinship), count in counts.items():
#         X.append(kinship)
#         Y.append(snp)
#         C.append(count)
#     X = np.array(X)
#     Y = np.array(Y)
#     C = np.array(C)
#     C2 = 75 * np.log2(C)
    
#     fig, ax = plt.subplots(1, figsize=(16,16))
#     ax.scatter(x=X,y=Y, s=C2, c='b', alpha=0.4)
#     ax.set_xlabel('Kinship')
#     ax.set_ylabel('# Variants')
#     fig.savefig(dest)


def plot_from_kinship_csv():
    kinship_results = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/kinship_results.csv'
    df = pd.read_csv(kinship_results, names=KINSHIP_CSV_HEADER, 
        sep=' ', usecols=['kinship', 'snps'])
    # df['person1'], df['person2'], _ = df.pair.str.split('\t').str
    # del df['pair']

    fig, ax = plt.subplots(figsize=(10,12))
    sns.distplot(df.kinship, kde=False, ax=ax)
    plt.yscale('log', nonposy='clip');
    plt.title('Histogram of Kinship')
    plt.savefig('/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/histogram_of_kinship.png')

    fig, ax = plt.subplots(figsize=(10,12))
    sns.distplot(df.snps, kde=False, ax=ax)
    plt.yscale('log', nonposy='clip');
    plt.title('Histogram of # SNPs')
    plt.savefig('/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/histogram_of_n_snps.png')

    sns.jointplot('snps', 'kinship', kind='hex', data=df, size=10)
    plt.savefig('/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/hex.png')
    # sns.lmplot('snps', 'kinship', data=df, size=6)
    # plt.savefig('lmplot.png')


    sns.jointplot('snps', 'kinship', data=df, s=10, size=10).plot_joint(sns.kdeplot, zorder=0, n_levels=6)
    plt.savefig('/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/kde.png')


def generate_missingness_stats(vcf, plink_bin, out):
    command = f'{plink_bin} --missing --vcf {vcf} --out {out}'
    run(command)
    path = f'{out}.lmiss'
    df = pd.read_csv(path, sep=r'\s+', engine='python')
    fig, ax = plt.subplots(1, figsize=(10,10))
    sns.distplot(df['F_MISS'], kde=False, bins=15, ax=ax)
    plt.xticks(np.arange(0, 1, 0.1))
    plt.xlabel('Missingness')
    plt.ylabel('# Variants')
    path = f'{out}.lmiss.png'
    plt.savefig(path)


def generate_LD_stats(vcf, plink_bin, out):
    command = f'{plink_bin} --r2 dprime with-freqs --vcf {vcf} --out {out}'
    run(command)


def generate_MAF_stats(vcf, plink_bin, out):
    n_samples = len(get_sample_ids_from_vcf(vcf))
    command = f'{plink_bin} --freq --vcf {vcf} --out {out}'
    run(command)
    path = f'{out}.frq'
    df = pd.read_csv(path, sep=r'\s+', engine='python')
    df['# Observed Alleles / (2 x # Samples)'] = df['NCHROBS'] / (2 * n_samples)
    sns.set_context('talk')
    g = sns.jointplot(x='MAF', y='# Observed Alleles / (2 x # Samples)', data=df, kind='kde', stat_func=None, xlim=(0, 0.5), ylim=(0, 1), size=12)
    g.plot_joint(plt.scatter, c='orange', s=40)
    g.ax_joint.collections[0].set_alpha(0)
    path = f'{out}.frq.png'
    g.savefig(path)
