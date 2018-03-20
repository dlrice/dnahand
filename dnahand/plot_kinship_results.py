#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white", color_codes=True)

def main():
    names = ['pair', 'ibd0', 'ibd1', 'ibd2', 'kinship', 'snps']
    kinship_results = '/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/dnahand-out/kinship/kinship_results.csv'
    df = pd.read_csv(kinship_results, names=names, sep=' ', usecols=['kinship', 'snps'])
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

if __name__ == '__main__':
    main()