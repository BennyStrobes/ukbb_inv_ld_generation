#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-15:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)





chrom_num="$1"
trait_summary_file="$2"
overlapping_windows_file="$3"
ukbb_in_sample_ld_dir="$4"
ukbb_in_sample_genotype_dir="$5"
processed_ukbb_in_sample_ld_dir="$6"


python3 generate_in_sample_ld_matrices.py $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir
