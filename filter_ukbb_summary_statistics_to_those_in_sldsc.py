import numpy as np
import os
import sys
import pdb
import gzip

def extract_dictionary_of_ukbb_pqtl_snps(ukbb_pqtl_genotype_data_dir):
	dicti = {}
	for chrom_num in range(1,23):
		chrom_dicti = {}
		chrom_pvar_file = ukbb_pqtl_genotype_data_dir + 'chr' + str(chrom_num) + '.pvar'
		f = open(chrom_pvar_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[2]
			if rsid.startswith('rs') == False and rsid.startswith('ss') == False:
				print('error: snp skipped cause didnt start with rs (ok if not too many of these errors)')
				continue
			if rsid in chrom_dicti:
				print('repeat rsid assumption erorro')
				pdb.set_trace()

			chrom_dicti[rsid] = (data[1], data[3], data[4])
		f.close()
		dicti[str(chrom_num)] = chrom_dicti
	return dicti


def extract_dictionary_of_ldsc_snps(ldsc_1kg_plink_dir):
	dicti = {}
	for chrom_num in range(1,23):
		chrom_dicti = {}
		chrom_plink_bim_file = ldsc_1kg_plink_dir + '1000G.EUR.QC.' + str(chrom_num) + '.bim'
		f = open(chrom_plink_bim_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			rsid = data[1]
			if rsid.startswith('rs') == False and rsid.startswith('ss') == False:
				print('error: snp skipped cause didnt start with rs (ok if not too many of these errors)')
				continue
			if rsid in chrom_dicti:
				print('repeat rsid assumption erorro')
				pdb.set_trace()

			chrom_dicti[rsid] = (data[3], data[4], data[5])
		f.close()
		dicti[str(chrom_num)] = chrom_dicti
	return dicti

def create_ordered_list_of_gwas_studies_by_parsing_input_directory(orig_ukbb_sumstats_hg19_dir):
	gwas_studies = []
	for file_name in os.listdir(orig_ukbb_sumstats_hg19_dir):
		if file_name.endswith('.bgen.stats.gz') == False:
			continue
		if file_name.endswith('.interim.bgen.stats.gz'):
			continue
		study_name = file_name.split('v3.')[-1].split('.bgen')[0]
		if study_name.startswith('460K'):
			continue
		gwas_studies.append(study_name)
	if len(gwas_studies) != len(np.unique(gwas_studies)):
		print('assumption eroror')
		pdb.set_trace()

	gwas_studies = np.sort(np.asarray(gwas_studies))
	return gwas_studies

def create_new_association_file(orig_gwas_file, new_gwas_file, ldsc_snp_dicti, ukbb_pqtl_snp_dicti):
	f = gzip.open(orig_gwas_file)
	t = open(new_gwas_file,'w')
	miss1 = 0
	miss2 = 0
	miss3 = 0
	head_count = 0
	counter = 0
	snps = {}
	missed_snps = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data[:8])) + '\t' + 'z_bolt_lmm\n')
			continue
		if len(data) != 16:
			print('assumption eroror')
			pdb.set_trace()
		chrom_num = data[1]
		rs_id = data[0]
		pos = data[2]
		if rs_id not in ldsc_snp_dicti[chrom_num] or rs_id not in ukbb_pqtl_snp_dicti[chrom_num]:
			continue
		counter = counter + 1
		a1 = data[4]
		a2 = data[5]
		(ldsc_snp_pos, ldsc_snp_a1, ldsc_snp_a2) = ldsc_snp_dicti[chrom_num][rs_id]
		(ukbb_snp_pos, ukbb_snp_a1, ukbb_snp_a2) = ukbb_pqtl_snp_dicti[chrom_num][rs_id]


		# Error checking
		if ldsc_snp_pos != pos:
			miss1 = miss1 + 1
			continue

		if ukbb_snp_pos != pos:
			print('assumption eroror')
			pdb.set_trace()

		if ldsc_snp_a1 != a1 or ldsc_snp_a2 != a2:
			if ldsc_snp_a1 != a2 or ldsc_snp_a2 != a1:
				missed_snps[rs_id] = 1
				continue
		if ukbb_snp_a1 != a1 or ukbb_snp_a2 != a2:
			if ukbb_snp_a1 != a2 or ukbb_snp_a2 != a1:
				miss3 = miss3 + 1
				continue


		if rs_id in snps:
			print('assumption eroror')
			pdb.set_trace()
		snps[rs_id] = 1

		beta = float(data[10])
		std_err = float(data[11])
		chi_sq_bolt_lmm = float(data[14])
		if chi_sq_bolt_lmm < 0.0:
			miss2 = miss2 + 1
			continue
		z_bolt_lmm = np.sqrt(chi_sq_bolt_lmm)
		signer = 1.0
		if beta <= 0.0:
			signer = -1.0
		z_bolt_lmm = z_bolt_lmm*signer
		t.write('\t'.join(np.asarray(data[:8])) + '\t' + '\t' + str(z_bolt_lmm) + '\n')
	t.close()
	return

def extract_list_of_snps_present_in_all_studies(ordered_gwas_studies, ukbb_sumstats_hg19_dir):
	dicti = {}
	for gwas_study in ordered_gwas_studies:
		temp_new_gwas_file = ukbb_sumstats_hg19_dir + gwas_study + '_hg19_sumstats_sldsc_snp_filtered_tmp.txt'
		f = open(temp_new_gwas_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			rsid = data[0]
			if rsid not in dicti:
				dicti[rsid] = (1, data[2], data[4], data[5])
			else:
				old_tupler = dicti[rsid]
				if old_tupler[1] != data[2] or old_tupler[2] != data[4] or old_tupler[3] != data[5]:
					print('assumption eroror')
					pdb.set_trace()
				dicti[rsid] = (old_tupler[0] + 1, data[2], data[4], data[5])
		f.close()

	final_dicti = {}
	for rsid in [*dicti]:
		if dicti[rsid][0] == len(ordered_gwas_studies):
			final_dicti[rsid] = 1
		elif dicti[rsid][0] > len(ordered_gwas_studies):
			print('rsid repeat error')
			pdb.set_trace()

	return final_dicti

def filter_association_file_to_new_snp_set(temp_new_gwas_file, new_gwas_file, snps_present_all_studies):
	f = open(temp_new_gwas_file)
	t = open(new_gwas_file,'w')
	head_count = 0
	used_snps = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		rsid = data[0]
		if rsid not in snps_present_all_studies:
			continue

		if rsid in used_snps:
			print('repeat snp assumption eroror')
			pdb.set_trace()
		used_snps[rsid] = 1
		t.write(line + '\n')
	f.close()
	t.close()
	return


#######################
# Command line args
#######################
orig_ukbb_sumstats_hg19_dir = sys.argv[1]
ukbb_sumstats_hg19_dir = sys.argv[2]
ldsc_1kg_plink_dir = sys.argv[3]
ukbb_pqtl_genotype_data_dir = sys.argv[4]


# Create ordered list of gwas studies
ordered_gwas_studies = create_ordered_list_of_gwas_studies_by_parsing_input_directory(orig_ukbb_sumstats_hg19_dir)
######*****##########


# Extract dictionary of snps used by ldsc (1KG)
ldsc_snp_dicti = extract_dictionary_of_ldsc_snps(ldsc_1kg_plink_dir)

# Extract dictionary list of snps used by UKBB PQTL data
ukbb_pqtl_snp_dicti = extract_dictionary_of_ukbb_pqtl_snps(ukbb_pqtl_genotype_data_dir)


# Create study files
study_organizer_file = ukbb_sumstats_hg19_dir + 'ukbb_hg19_sumstat_files.txt'
t_outer = open(study_organizer_file,'w')
t_outer.write('study_name\tstudy_file\n')
# Loop through studies

for itera, gwas_study in enumerate(ordered_gwas_studies):
	print(gwas_study)
	# Original gwas file in hg19 coordinate space
	orig_gwas_file = orig_ukbb_sumstats_hg19_dir + 'bolt_337K_unrelStringentBrit_MAF0.001_v3.' + gwas_study + '.bgen.stats.gz'

	# First convert hg19 association to temporary bed format
	temp_new_gwas_file = ukbb_sumstats_hg19_dir + gwas_study + '_hg19_sumstats_sldsc_snp_filtered_tmp.txt'
	new_gwas_file = ukbb_sumstats_hg19_dir + gwas_study + '_hg19_sumstats_sldsc_snp_filtered.txt'
	create_new_association_file(orig_gwas_file, temp_new_gwas_file, ldsc_snp_dicti, ukbb_pqtl_snp_dicti)


	# Print study and new lifted over file to output
	t_outer.write(gwas_study + '\t' + new_gwas_file + '\n')

t_outer.close()



snps_present_all_studies = extract_list_of_snps_present_in_all_studies(ordered_gwas_studies, ukbb_sumstats_hg19_dir)

for itera, gwas_study in enumerate(ordered_gwas_studies):
	print(gwas_study)
	# First convert hg19 association to temporary bed format
	temp_new_gwas_file = ukbb_sumstats_hg19_dir + gwas_study + '_hg19_sumstats_sldsc_snp_filtered_tmp.txt'
	new_gwas_file = ukbb_sumstats_hg19_dir + gwas_study + '_hg19_sumstats_sldsc_snp_filtered.txt'
	filter_association_file_to_new_snp_set(temp_new_gwas_file, new_gwas_file, snps_present_all_studies)
	
	# Remove temp file
	os.system('rm ' + temp_new_gwas_file)
