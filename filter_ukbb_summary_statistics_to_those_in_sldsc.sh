#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=15G                         # Memory total in MiB (for all cores)


orig_ukbb_sumstats_hg19_dir="$1"
ukbb_sumstats_hg19_dir="$2"
ldsc_1kg_plink_dir="$3"
ukbb_pqtl_genotype_data_dir="$4"

date
python3 filter_ukbb_summary_statistics_to_those_in_sldsc.py ${orig_ukbb_sumstats_hg19_dir} ${ukbb_sumstats_hg19_dir} ${ldsc_1kg_plink_dir} ${ukbb_pqtl_genotype_data_dir}
date

python3 generate_trait_list_file_with_sample_size_and_heritabilities.py ${ukbb_sumstats_hg19_dir} ${orig_ukbb_sumstats_hg19_dir}
