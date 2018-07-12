from download import download_fingerprints, read_sangerids, digest_sample_lists_directory
from handprint import fluidigm, sequenom
from analysis import generate_missingness_stats, generate_LD_stats, generate_MAF_stats, get_counts_from_csv
import utils
import os
import sys
from glob import glob
import pickle

def download_collate_to_vcf_kinship(
    sample_lists_to_download_directory, out_directory, 
    reference_vcf_path, reference_snp_pickle, chromosomes, 
    vcf_from_plex_bin, bcftools_bin, baton_bin, baton_metaquery_bin, 
    baton_get_bin, akt, info_include_path, irods_credentials_path=None,
    n_max_processes=25, plink_bin=None, pipeline_entry_name='download'):
    """
    Args:
        sample_list_path (str): Path to a headerless text file which
            lists the Sanger sample IDs to process.
        out_directory (str): Directory where the all results will be
            saved. If already exists returns error.
        irods_credentials_path (str, optional): Path to a text file
            containing user's irods password. If not supplied, user will
            be prompted.
        reference_vcf_path (str): Path to VCF which will be merged
            against.
        vcf_from_plex_bin (str):
        bcftools_bin (str):
        entry_point (str): The point at which the pipeline starts:
            -download
            -makehandprints
            -makevcf
            -kinship

    Returns:
        Creates the following directory structure upon completion:
            out_directory/
                execution_arguments.txt
                fingerprints/sequenom: downloaded CSV files
                fingerprints/fluidigm: downloaded CSV files
                handprints/sequenom/handprint*/
                    sampleid*.csv
                    snpset.tsv
                handprints/fluidigm/handprint*/
                    sampleid*.csv
                    snpset.tsv
                logs/vcf_from_plex_fluidigm_handprint1.o|e
                vcfs/
                    sequenom/
                        handprint*.vcf.gz
                    fluidigm/
                        handprint*.vcf.gz
                    handprints.vcf.gz(.csi)
                    subsetted_reference.vcf.gz(.csi)
                    all_merged.vcf.gz(.csi)
                kinship/kinship_results.csv
                kinship/kinship_results.pickle
    """


    """
    TODO
    1. Determine how much the handprints stuff actually helps. How many
       times does one individual have multiple fingerprints?
    2. Plot 
    3. Prune on LD in reference vcf
    """

    utils.PLINK_BIN = plink_bin
    utils.BCFTOOLS_BIN = bcftools_bin

    PIPELINE_STEPS = {
        'download': 0,
        'makehandprints': 1,
        'plex2vcf': 2,
        'mergevcfs': 3,
        'check_sampleids': 4,
        'kinship': 5,
        'analysis': 6,
    }

    assert pipeline_entry_name in PIPELINE_STEPS, \
        'Please specify a valid entry point.'

    pipeline_entry = PIPELINE_STEPS.get(pipeline_entry_name, 0)

    print(f'Starting at {pipeline_entry_name} in pipeline.')

    fingerprints_directory = os.path.join(out_directory, 'fingerprints')

    if pipeline_entry <= PIPELINE_STEPS['download']:
        # Make sure directory doesn't already exist - we want to start from
        # fresh.
        # if os.path.exists(out_directory):
        #     print(f'Error: {out_directory} already exists. Please provide '
        #         'a path to nonexistent directory to deposit results.')
        #     sys.exit(-1)

        # Download raw fingerprints from iRODS
        os.makedirs(out_directory, exist_ok=True)
        os.makedirs(fingerprints_directory, exist_ok=True)
        
        sample_list_irods_db_path = os.path.join(out_directory, 'sample_list_irods.json')

        digest_sample_lists_directory(
            sample_lists_to_download_directory, sample_list_irods_db_path,
            fingerprints_directory, baton_bin,
            baton_metaquery_bin, baton_get_bin, irods_credentials_path,
            n_max_processes)
        # for fingerprint_method in utils.FINGERPRINT_METHODS:
        #     # Create directory for method
        #     fingerprint_method_directory = os.path.join(
        #         fingerprint_directory,fingerprint_method)
        #     os.makedirs((fingerprint_method_directory, exist_ok=True)
        #     download_fingerprints(sample_list_path, 
        #         fingerprints_directory, fingerprint_method,
        #         baton_bin, baton_metaquery_bin, baton_get_bin,
        #         n_max_processes)

    handprint_directory = os.path.join(out_directory, 'handprints')
    # Generate handprints from downloaded
    if pipeline_entry <= PIPELINE_STEPS['makehandprints']:
        os.mkdir(handprint_directory)
        for fingerprint_method in utils.FINGERPRINT_METHODS:
            fingerprint_method_directory = os.path.join(
                fingerprint_directory,
                fingerprint_method)
            handprint_method_directory = os.path.join(
                handprint_directory,
                fingerprint_method)
            if fingerprint_method == 'fluidigm':
                fluidigm.generate(fingerprint_method_directory,
                    handprint_method_directory, reference_snp_pickle)
            elif fingerprint_method == 'sequenom':
                sequenom.generate(fingerprint_method_directory,
                    handprint_method_directory, reference_snp_pickle)

    vcf_directory = os.path.join(out_directory, 'vcfs')
    print(vcf_directory)
    # Convert handprints to VCFs
    if pipeline_entry <= PIPELINE_STEPS['plex2vcf']:
        os.mkdir(vcf_directory)
        # log_directory = os.path.join(out_directory, 'logs')
        # os.mkdir(log_directory)
        for fingerprint_method in utils.FINGERPRINT_METHODS:
            handprint_method_directory = os.path.join(handprint_directory,
                fingerprint_method)
            handprints = utils.get_subdirectories(
                handprint_method_directory)

            n_handprints = len(handprints)
            if not n_handprints:
                print(f'No handprints generated for {fingerprint_method}')
                continue

            vcf_out_dir = os.path.join(vcf_directory, fingerprint_method)
            if not os.path.exists(vcf_out_dir):
                os.makedirs(vcf_out_dir)

            for handprint in handprints:
                filelist = os.path.join(handprint, f'filelist.txt')
                
                snpset = os.path.join(handprint, 'snpset.tsv')
                handprint = os.path.normpath(handprint)
                handprint = os.path.basename(handprint)
                vcf_filename = f'{handprint}.vcf'
                vcf_out = os.path.join(vcf_out_dir, vcf_filename)
                print(vcf_out)

                # # Save output of vcf_from_plex to log files
                # _, this_handprint_directory = os.path.split(handprint)
                # base = f'{fingerprint_method}_{this_handprint_directory}'
                # stdout = os.path.join(log_directory, base + '.o')
                # stderr = os.path.join(log_directory, base + '.e')

                utils.run_vcf_from_plex(vcf_from_plex_bin, filelist,
                    chromosomes, fingerprint_method, snpset, vcf_out)
                utils.convert_sampleids_to_lowercase_vcf(vcf_out)

    # Merge all VCFs
    vcf_merged_path = os.path.join(vcf_directory, 'all_merged.vcf.gz')
    subsetted_reference_vcf = os.path.join(vcf_directory, 
        'subsetted_reference.vcf.gz')
    vcf_merged_handprints_path = os.path.join(vcf_directory, 
        'handprints.vcf.gz')
    if pipeline_entry <= PIPELINE_STEPS['mergevcfs']:
        handprint_vcfs = glob(os.path.join(vcf_directory, '*', '*'))
        print(handprint_vcfs)
        n_handprint_vcfs = len(handprint_vcfs)
        if not n_handprint_vcfs:
            print('Error: No handprint VCFs found to merge.')
            sys.exit(-1)

        print(f'Found {n_handprint_vcfs} handprint VCFs to merge')

        vcf_subset_path = os.path.join(vcf_directory,
        f'subsetted_reference.vcf.gz')
        subsetted_reference_vcf = os.path.join(vcf_directory, 
            'subsetted_reference.vcf.gz')
        # Subset reference on only those SNPs seen in all of the 
        # handprint_vcfs
        utils.subset_reference_vcfs_on_handprint_snps(handprint_vcfs,
            reference_snp_pickle, reference_vcf_path,
            subsetted_reference_vcf)
        utils.index_vcf(subsetted_reference_vcf)
        
        # utils.concat_vcfs(subsetted_referenence_vcfs, referenence_vcf)
        
        
        # Remove samples in handprints already present in reference - these
        # will be much higher quality in the imputed reference
        handprint_vcfs = utils.filter_duplicate_individuals(
            subsetted_reference_vcf, handprint_vcfs)
        handprint_vcfs = utils.gzip_vcfs(handprint_vcfs)
        utils.index_vcfs(handprint_vcfs)

        # Create a VCF of all handprints to print MAF/LD/missing stats
        utils.merge_vcfs(handprint_vcfs, vcf_merged_handprints_path)
        utils.index_vcf(vcf_merged_handprints_path)

        vcfs_to_merge = [
            subsetted_reference_vcf,
            vcf_merged_handprints_path
        ]

        # Merge all VCFs
        utils.merge_vcfs(vcfs_to_merge, vcf_merged_path)
        utils.index_vcf(vcf_merged_path)

        # Prune subsetted_reference_vcf by LD - don't trust the
        # fingerprints because of the (most likely) high missingness.
        prune_exclude_path = utils.get_LD_prune_list(
            subsetted_reference_vcf)
        print('Pruning the following:')
        with open(prune_exclude_path) as f:
            print(f.readlines())
        utils.prune_by_rsid(vcf_merged_path, prune_exclude_path)
        # Post-pruning will still have the same filename
        utils.index_vcf(vcf_merged_path)

        # Generate stats
        vcfs = vcfs_to_merge + [vcf_merged_path]

        for vcf in vcfs:
            i = vcf.rfind('.vcf')
            out = vcf[:i]
            generate_missingness_stats(vcf, plink_bin, out)
            generate_LD_stats(vcf, plink_bin, out)
            generate_MAF_stats(vcf, plink_bin, out)

        # Filter out low MAF variants
        i = vcf_merged_path.rfind('.vcf')
        root = vcf_merged_path[:i] 
        frq_path = f'{root}.frq'
        maf_out = f'{root}.maf.out' 
        utils.get_MAF_exclude_list(frq_path, maf_out)
        utils.prune_by_rsid(vcf_merged_path, maf_out)
        utils.index_vcf(vcf_merged_path)


        # Filter out low INFO variants
        utils.include_by_region_file(vcf_merged_path, info_include_path)
        utils.index_vcf(vcf_merged_path)
    
    if pipeline_entry <= PIPELINE_STEPS['check_sampleids']:
        print('Checking number of the Sanger sample IDs present in the final merged VCF.')
        sample_ids_for_lookup = read_sangerids(sample_list_path)
        sample_ids_for_lookup = {sampleid.lower() for sampleid in sample_ids_for_lookup if sampleid}
        n_sample_ids_for_lookup = len(sample_ids_for_lookup)
        sample_ids_from_vcf = utils.get_sample_ids_from_vcf(vcf_merged_path)
        n_included = len(sample_ids_from_vcf & sample_ids_for_lookup)
        n_not_included = len(sample_ids_for_lookup) - n_included
        print(f'...N included {n_included} / {n_sample_ids_for_lookup}')
        print(f'...N not included {n_not_included} / {n_sample_ids_for_lookup}')
        assert n_included == n_sample_ids_for_lookup

    # Run akt
    kinship_directory = os.path.join(out_directory, 'kinship')
    kinship_results = os.path.join(kinship_directory, 'kinship_results_minkin_0.1.csv')
    if pipeline_entry <= PIPELINE_STEPS['kinship']:
        os.makedirs(kinship_directory, exist_ok=True)
        kinship_results = os.path.join(kinship_directory, 'kinship_results_minkin_0.4.csv')
        command = f'{akt} kin --method 0 {vcf_merged_path} --force --minkin {0.4} > {kinship_results}'
        utils.run(command)

    analysis_directory = os.path.join(out_directory, 'analysis')
    if pipeline_entry <= PIPELINE_STEPS['analysis']:
        os.makedirs(analysis_directory, exist_ok=True)
        reference_samples = utils.get_sample_ids_from_vcf(subsetted_reference_vcf)
        handprint_samples = utils.get_sample_ids_from_vcf(vcf_merged_handprints_path)
        get_counts_from_csv(kinship_results, reference_samples, handprint_samples, analysis_directory)
