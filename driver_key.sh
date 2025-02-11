##################
# Input data
##################

# Directory containing summary statistics
orig_ukbb_sumstats_hg19_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"

# Hapmap3 rsids
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"


# UKBB in sample LD
# Generated by Martin (though this file is temporary)
ukbb_in_sample_ld_dir="/n/groups/price/UKBiobank/insample_ld/ld_files/"
ukbb_in_sample_genotype_dir="/n/groups/price/UKBiobank/insample_ld/pvar_files/"

# LDSC baseline LD Dir
ldsc_baseline_ld_hg19_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD_v2.2/"

# 1KG snps used in SLDSC
ldsc_1kg_plink_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/"

# UKBB PQTL Genotype data
ukbb_pqtl_genotype_data_dir="/n/scratch/users/k/kah3107/UKB-PPP-DATA/genotype/imputed/"


##################
# Output data
##################
# Output root directory
output_root="/n/scratch/users/b/bes710/ukbb_inv_ld/"
perm_output_root="/n/groups/price/ben/ukbb_inv_ld/"

# Directory containing filtered summary statistics
ukbb_sumstats_hg19_dir=$output_root"ukbb_sumstats_hg19/"

# Directory containing list of 3MB windows
overlapping_3_mb_windows_dir=${perm_output_root}"overlapping_3mb_windows/"

# Directory containing in sample LD matrices
processed_ukbb_in_sample_ld_dir=${output_root}"processed_ukbb_in_sample_ld/"



########################################
########################################
# Analysis
########################################


# Filter UKBB summary statistics to those in SLDSC (1KG)
if false; then
sbatch filter_ukbb_summary_statistics_to_those_in_sldsc.sh ${orig_ukbb_sumstats_hg19_dir} ${ukbb_sumstats_hg19_dir} ${ldsc_1kg_plink_dir} ${ukbb_pqtl_genotype_data_dir}
fi
trait_summary_file=${ukbb_sumstats_hg19_dir}"ukbb_hg19_sumstat_files_with_samp_size_and_h2.txt"


# Generate overlapping 3MB LD windows
overlapping_windows_file=${overlapping_3_mb_windows_dir}"genome_wide_windows.txt"
if false; then
sh generate_overlapping_3mb_ld_windows.sh $ukbb_sumstats_hg19_dir $overlapping_windows_file
fi


# Generate in sample LD in each chromosome (do seperately for each chromosome)
# Note that chrom 8 needed 45GB of memory (where rest need 20GB)
if false; then
for chrom_num in $(seq 1 22); do 
	sbatch generate_in_sample_ld_matrices.sh $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir
done
fi

if false; then
chrom_num="1"
sbatch generate_in_sample_ld_matrices.sh $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir

chrom_num="8"
sbatch generate_in_sample_ld_matrices.sh $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir

chrom_num="9"
sbatch generate_in_sample_ld_matrices.sh $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir

chrom_num="16"
sbatch generate_in_sample_ld_matrices.sh $chrom_num $trait_summary_file $overlapping_windows_file $ukbb_in_sample_ld_dir $ukbb_in_sample_genotype_dir $processed_ukbb_in_sample_ld_dir
fi

