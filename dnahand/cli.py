#!/usr/bin/env python3
from pipeline import download_collate_to_vcf_kinship
import handprint.snp_reference
import argparse

def main():
    # download_collate_to_vcf_kinship(
    #     sample_list_path='',
    #     out_directory='/lustre/scratch115/teams/barrett/users/dr9/dnahand-out-complete',
    #     reference_vcf_path='/lustre/scratch115/realdata/mdt0/teams/barrett/users/dr9/new_imputation_dan/all_studies.sampleids_cleaned_to_lowercase.bcf.gz',
    #     reference_snp_pickle='/nfs/team143/fingerprint_resources/all_rsids_in_imputed_datasets.pickle',
    #     chromosomes='/nfs/team143/fingerprint_resources/chromosome_lengths_GRCh37.json',
    #     vcf_from_plex_bin='module load genotyping/1.14.2; /software/gapi/pkg/genotyping/1.14.2/bin/vcf_from_plex.pl',
    #     bcftools_bin='/software/hgi/pkglocal/bcftools-1.6-htslib-1.6-htslib-plugins-6f2229e0-irods-git-4.2.2-plugin_kerberos-2.0.0/bin-wrap/bcftools',
    #     baton_bin='/software/solexa/pkg/baton/0.17.1/bin/baton',
    #     baton_metaquery_bin='/software/solexa/pkg/baton/0.17.1/bin/baton-metaquery',
    #     baton_get_bin='/software/solexa/pkg/baton/0.17.1/bin/baton-get',
    #     irods_credentials_path='/nfs/users/nfs_d/dr9/.temp',
    #     pipeline_entry_name='mergevcfs',
    #     akt='/nfs/team143/akt/akt'
    # )

    parser = argparse.ArgumentParser()
    parser.add_argument('task',
        choices=['pipeline', 'save_kinship_frequency_pickle', 'generate_snp_data_pickle'])

    # pipeline arguments
    parser.add_argument('--sample_list_path')
    parser.add_argument('--out_directory')
    parser.add_argument('--reference_vcf_path')
    parser.add_argument('--reference_snp_pickle')
    parser.add_argument('--chromosomes')
    parser.add_argument('--vcf_from_plex_bin')
    parser.add_argument('--bcftools_bin')
    parser.add_argument('--baton_bin')
    parser.add_argument('--baton_metaquery_bin')
    parser.add_argument('--baton_get_bin')
    parser.add_argument('--irods_credentials_path')
    parser.add_argument('--pipeline_entry_name')
    parser.add_argument('--akt')
    parser.add_argument('--plink')

    # generate_snp_data_pickle arguments
    # --reference-vcf-path
    parser.add_argument('--snp_data_dest')
    parser.add_argument('--snp_data_pickle_dest')

    args = parser.parse_args()


    if args.task == 'generate_snp_data_pickle':
        # handprint.snp_reference.create_snp_data_textfile(
        #     reference_vcf_path=args.reference_vcf_path,
        #     snp_data_dest=args.snp_data_dest, 
        #     bcftools_bin=args.bcftools_bin
        # )
        snp_data = handprint.snp_reference.read_snp_data_from_textfile(
            snp_data_path=args.snp_data_dest)
        handprint.snp_reference.save_snp_data_as_pickle(
            snp_data=snp_data,
            snp_data_pickle_dest=args.snp_data_pickle_dest
        )
    elif args.task == 'pipeline':
        download_collate_to_vcf_kinship(
            sample_list_path=args.sample_list_path,
            out_directory=args.out_directory,
            reference_vcf_path=args.reference_vcf_path,
            reference_snp_pickle=args.reference_snp_pickle,
            chromosomes=args.chromosomes,
            vcf_from_plex_bin=args.vcf_from_plex_bin,
            bcftools_bin=args.bcftools_bin,
            baton_bin=args.baton_bin,
            baton_metaquery_bin=args.baton_metaquery_bin,
            baton_get_bin=args.baton_get_bin,
            irods_credentials_path=args.irods_credentials_path,
            pipeline_entry_name=args.pipeline_entry_name,
            akt=args.akt,
            plink=args.plink
        )


if __name__ == '__main__':
    main()