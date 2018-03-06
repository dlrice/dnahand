from download import 
from handprint import handprint
import utils
import os
import sys

def download_collate_to_vcf_kinship(sample_list_path, out_directory, 
    reference_vcf_pattern, reference_snp_pickle, chromosomes, 
    vcf_from_plex_bin, bcftools_binbcftools_bin, baton_bin, baton_metaquery_bin, 
    baton_get_bin, irods_credentials_path=None, n_max_processes=25,
    entry_point=None):
    """
    Args:
        sample_list_path (str): Path to a headerless text file which
            lists the Sanger sample IDs to process.
        out_directory (str): Directory where the all results will be
            saved. If already exists returns error.
        irods_credentials_path (str, optional): Path to a text file
            containing user's irods password. If not supplied, user will
            be prompted.
        reference_vcf_pattern (str): Path pattern which when supplied
            with a chromosome: reference_vcf_pattern.format(chromosome)
            points to a VCF which will be merged against
        vcf_from_plex_bin (str):
        bcftools_bin (str):
        entry_point (str): The point at which the pipeline starts:
            -download
            -makehandprints
            -makevcf
            -runkinship

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
                    reference/
                        reference*.vcf.gz
                        reference.vcf.gz
                    all_merged.vcf.gz

                kinship/kinship_results.csv
                kinship/kinship_results.pickle
    """

    utils.BCFTOOLS_BIN = bcftools_bin
    PIPELINE_STEPS = [
        'download': 0,
        'makehandprints': 1,
        'plex2vcf': 2,
        'mergevcfs': 3
        'runkinship': 4,
    ]

    assert (PIPELINE_STEPS is None) or (entry_point in PIPELINE_STEPS),
        'Please specify a valid entry point.'

    pipeline_start = PIPELINE_STEPS.get(entry_point, 0)

    fingerprint_directory = os.path.join(out_directory, 'fingerprints')

    if pipeline_start <= PIPELINE_STEPS['download']:
        # Make sure directory doesn't already exist - we want to start from
        # fresh.
        if os.path.exists(out_directory):
            print(f'Error: {out_directory} already exists. Please provide '
                'a path to nonexistent directory to deposit results.')
            sys.exit(-1)

        # Download raw fingerprints from iRODS
        os.mkdir(out_directory)
        os.mkdir(fingerprint_directory)
        for fingerprint_method in FINGERPRINT_METHODS:
            # Create directory for method
            fingerprint_method_directory = os.path.join(
                fingerprint_directory,fingerprint_method)
            os.mkdir(fingerprint_method_directory)
            download_fingerprints(sample_list_path,
                fingerprint_method_directory, fingerprint_method, 
                n_max_processes, baton_bin, baton_metaquery_bin, 
                baton_get_bin)

    # Generate handprints from downloaded
    if pipeline_start <= PIPELINE_STEPS['makehandprints']:
        handprint_directory = os.path.join(out_directory, 'handprints')
        os.mkdir(handprint_directory)
        for fingerprint_method in FINGERPRINT_METHODS:
            fingerprint_method_directory = os.path.join(
                fingerprint_directory,
                fingerprint_method)
            handprint_method_directory = os.path.join(handprint_directory,
                fingerprint_method)
            if fingerprint_method == 'fluidigm':
                fluidigm.generate(fingerprint_method_directory,
                    handprint_method_directory, reference_snp_pickle)
            elif fingerprint_method == 'sequenom':
                sequenom.generate(fingerprint_method_directory,
                    handprint_method_directory, reference_snp_pickle)


    vcf_directory = os.path.join(out_directory, 'vcfs')
    # Convert handprints to VCFs
    if pipeline_start <= PIPELINE_STEPS['plex2vcf']:
        os.mkdir(vcf_directory)
        log_directory = os.path.join(out_directory, 'logs')
        os.mkdir(log_directory)
        for fingerprint_method in FINGERPRINT_METHODS:
            handprint_method_directory = os.path.join(handprint_directory,
                fingerprint_method)
            handprints = utils.get_subdirectories(
                handprint_method_directory)
            n_handprints = len(handprints)
            if not n_handprints:
                print(f'No handprints generated for {fingerprint_method}')
                continue
            for handprint in handprints:
                filelist = os.path.join(handprint, f'filelist.txt')
                
                snpset = os.path.join(handprint, 'snpset.tsv')
                vcf_filename = f'{handprint}.vcf.gz'
                vcf_out = os.path.join(vcf_directory, fingerprint_method,
                    vcf_filename)

                # Save output of vcf_from_plex to log files
                _, this_handprint_directory = os.path.split(handprint)
                base = f'{fingerprint_method}_{this_handprint_directory}'
                stdout = os.path.join(log_directory, base + '.o')
                stderr = os.path.join(log_directory, base + '.e')

                utils.run_vcf_from_plex(vcf_from_plex_bin, filelist,
                    chromosomes, fingerprint_method, snpset, vcf_out)

    # Merge all VCFs
    if pipeline_start <= PIPELINE_STEPS['mergevcfs']:
        handprint_vcfs = glob(vcf_directory, '*', '*')
        n_handprint_vcfs = len(handprint_vcfs)
        if not n_handprint_vcfs:
            print('Error: No handprint VCFs found to merge.')
            sys.exit(-1)

        reference_vcf_subset_directory = os.path.join(vcf_directory, 
            'reference')
        # Subset reference on only those SNPs seen in all of the 
        # handprint_vcfs
        utils.subset_reference_vcfs_on_handprint_snps(handprint_vcfs, 
            reference_snp_pickle, reference_vcf_pattern, 
            reference_vcf_subset_directory)
        subsetted_referenence_vcfs = glob(reference_vcf_subset_directory,
            '*')
        n_reference_vcfs = len(reference_vcfs)
        print('Found {n_handprint_vcfs} handprint VCFs to merge')
        print('Found {n_reference_vcfs} reference VCFs to merge')
        referenence_vcf = os.path.join(vcf_directory, 'reference.vcf.gz')
        utils.concat_vcfs(subsetted_referenence_vcfs, referenence_vcf)
        vcfs_to_merge = [referenence_vcf] + handprint_vcfs
        
        # Remove samples in handprints already present in reference - these
        # will be much higher quality in the imputed reference
        vcfs_to_merge = utils.sort_by_number_of_snps(vcfs_to_merge)
        vcfs_to_merge = utils.filter_duplicate_individuals(vcfs_to_merge)

        # Merge all VCFs
        vcf_merged_path = os.path.join(vcf_directory, 'all_merged.vcf.gz')
        utils.merge_vcfs(vcfs_to_merge, vcf_merged_path)

    # Run akt



def main():
    download.print_message('hi')
    handprint.handprint()

if __name__ == '__main__':
    main()
