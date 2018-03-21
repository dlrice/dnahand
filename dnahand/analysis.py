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

KINSHIP_CSV_HEADER = ['pair', 'ibd0', 'ibd1', 'ibd2', 'kinship', 'snps']


def get_frequencies_from_csv(kinship_results_csv_path, pickle_dest):
    counts = defaultdict(int)
    with open(kinship_results_csv_path) as f:
        reader = DictReader(f, delimiter=' ',
            fieldnames=KINSHIP_CSV_HEADER)
        for row in reader:
            kinship = row['kinship'][:4]
            snps = row['snps']
            counts[(snps, kinship)]+=1

    counts_numeric_keys = defaultdict(int)
    for k, v in counts.items():
        snps = int(k[0])
        kinship = round(float(k[1]), 2)
        k = (snps, kinship)
        counts_numeric_keys[k] += v

    return counts_numeric_keys


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


def generate_missingness_stats(vcf, plink, out):
    command = f'{plink} --missing --vcf {vcf} --out {out}'
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


def generate_LD_stats(vcf, plink, out):
    command = f'{plink} --r2 dprime with-freqs --vcf {vcf} --out {out}'
    run(command)


def generate_MAF_stats(vcf, plink, out):
    n_samples = len(get_sample_ids_from_vcf(vcf))
    command = f'{plink} --freq --vcf {vcf} --out {out}'
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