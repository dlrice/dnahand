#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", color_codes=True)
from collections import defaultdict


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
        kinship = round(float(k[0]), 2)
        snps = int(k[1])
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

def main():

if __name__ == '__main__':
    main()