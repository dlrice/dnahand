/nfs/team152/dnahand/dnahand/cli.py \
update \
--sample_lists_to_download_directory /lustre/scratch119/humgen/teams/anderson/projects/dnahand-results/all/sample_lists_to_download/  \
--out_directory /lustre/scratch119/humgen/teams/anderson/projects/dnahand-results/all \
--reference_vcf_path /lustre/scratch115/realdata/mdt0/projects/ibdgwas/post_imputation/Dan_Merged_Studies/all_studies.sampleids_cleaned_to_lowercase.bcf.gz \
--reference_snp_pickle /lustre/scratch115/realdata/mdt0/projects/ibdgwas/post_imputation/Dan_Merged_Studies/all_studies.sampleids_cleaned_to_lowercase.snp_data.pickle \
--info_include_path /lustre/scratch115/realdata/mdt0/projects/ibdgwas/post_imputation/Dan_Merged_Studies/hrc_imputed_sites_info_0.9.tsv \
--chromosomes /nfs/team152/dnahand/chromosome_lengths_GRCh37.json \
--vcf_from_plex_bin "module load genotyping/1.14.2; /software/gapi/pkg/genotyping/1.14.2/bin/vcf_from_plex.pl" \
--bcftools_bin /software/hgi/pkglocal/bcftools-1.6-htslib-1.6-htslib-plugins-6f2229e0-irods-git-4.2.2-plugin_kerberos-2.0.0/bin-wrap/bcftools \
--baton_bin /software/solexa/pkg/baton/0.17.1/bin/baton \
--baton_metaquery_bin /software/solexa/pkg/baton/0.17.1/bin/baton-metaquery \
--baton_get_bin /software/solexa/pkg/baton/0.17.1/bin/baton-get \
--akt /nfs/team152/akt/akt \
--plink_bin /software/hgi/pkglocal/plink-1.90b4/bin/plink \
--pipeline_entry_name plex2vcf \

