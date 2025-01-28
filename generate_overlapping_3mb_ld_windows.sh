#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in


ukbb_sumstats_hg19_dir="${1}"
overlapping_windows_file="${2}"





python3 generate_overlapping_3mb_ld_windows.py $ukbb_sumstats_hg19_dir $overlapping_windows_file
