import sys
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import scipy.sparse
import time


def create_mapping_from_rsid_to_in_sample_variant_index(chrom_pvar_file):
	f = open(chrom_pvar_file)
	dicti = {}
	rs_id_to_alleles = {}
	rs_id_to_pos = {}
	head_count = 0
	indexer = 0
	ordered_rsids = []
	ordered_positions = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[2]
		alleles = data[4] + '_' + data[3]
		pos = data[1]
		if rsid in rs_id_to_alleles:
			print('assumption eroror')
			pdb.set_trace()
		rs_id_to_alleles[rsid] = (data[4], data[3])
		if rsid in rs_id_to_pos:
			print('assumption eroror')
			pdb.set_trace()
		rs_id_to_pos[rsid] = pos
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[rsid] = indexer
		indexer = indexer + 1
		ordered_rsids.append(rsid)
		ordered_positions.append(int(pos))
	f.close()
	return dicti, rs_id_to_alleles, rs_id_to_pos, np.asarray(ordered_rsids), np.asarray(ordered_positions)


def create_summary_statistic_rsid_dictionary(trait_summary_file, chrom_num):
	# Extract name of first trait
	# doesn't matter if its first or second trait
	f = open(trait_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		trait_file = data[1]
		break
	f.close()

	# Create mapping
	dicti = {}
	f = open(trait_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] != chrom_num:
			continue
		rsid = data[0]
		snp_pos = data[2]
		a1 = data[4]
		a2 = data[5]
		if rsid in dicti:
			print('assumption eroror')
			pdb.set_trace()
		dicti[rsid] = (snp_pos, a1, a2)
	f.close()
	return dicti

def indices_dont_lie_in_file(sample_ld_variant_indices, file_index_start, file_index_end):
	booler = True
	for indexer in sample_ld_variant_indices:
		if indexer >= file_index_start and indexer < file_index_end:
			booler = False
	return booler

def read_ld(fpath):
	"""
	Read LD files.
		- `_fullld.npy` : full_ld matrix, np.array(dtype=np.float32)
		- `_ld.npz` : ld matrix with SNPs in 10MB window, sp.sparse.csc_matrix(dtype=np.float32)
	Parameters
	----------
	fpath: str
		LD file path.
	Returns
	-------
	mat_ld : np.array(dtype=np.float32) or sp.sparse.csc_matrix(dtype=np.float32)
		LD matrix of dimension (n_ref_snp, n_snp)
	dic_range : dict
		
		- dic_range['chr'] : chromosome
		- dic_range['start'] : start position
		- dic_range['end'] : end position
		- dic_range['chr_ref'] : reference chromosome list (List)      
	"""
	
	# Check fpath
	err_msg = "fpath should end with one of ['_fullld.npy', '_ld.npz'] : %s" % fpath
	assert fpath.endswith("_fullld.npy") | fpath.endswith("_ld.npz"), err_msg
	
	if fpath.endswith("_fullld.npy"):
		mat_ld = np.load(fpath)
		temp_str = [x for x in fpath.split('.') if x.endswith('_fullld')][0]
		dic_range = parse_snp_range(temp_str)
		
	if fpath.endswith("_ld.npz"):
		mat_ld = scipy.sparse.load_npz(fpath)
		temp_str = [x for x in fpath.split('.') if x.endswith('_ld')][0]
		dic_range = parse_snp_range(temp_str)

	return mat_ld,dic_range

def parse_snp_range(snp_range):
	"""Get range of SNPs to analyze. 
	
	Parameters
	----------
	snp_range: str
		Example: 'c1_s0_e2000_r1'
	Returns
	-------
	dic_range : dict
		
		- dic_range['chr'] : chromosome
		- dic_range['start'] : start position
		- dic_range['end'] : end position
		- dic_range['chr_ref'] : reference chromosome list (List)
	"""

	dic_range = {x: None for x in ["chr", "start", "end", "chr_ref"]}

	for x in snp_range.split("_"):

		if x[0] == "c":
			dic_range["chr"] = int(x.replace("c", "").strip())

		if x[0] == "s":
			dic_range["start"] = int(x.replace("s", "").strip())

		if x[0] == "e":
			dic_range["end"] = int(x.replace("e", "").strip())

		if x[0] == "r":
			temp_str = x.replace("r", "").strip()
			if temp_str == "all":
				dic_range["chr_ref"] = list(np.arange(1, 23))
			else:
				dic_range["chr_ref"] = [int(x) for x in temp_str.split(",")]

	return dic_range


def extract_ld_mat_from_in_sample_ld(sample_ld_variant_indices, ukbb_in_sample_ld_dir, chrom_num):
	min_index = np.min(sample_ld_variant_indices)
	max_index = np.max(sample_ld_variant_indices)

	num_var = len(sample_ld_variant_indices)

	ld_mat = np.zeros((num_var, num_var)) + -2000.0

	for file_name in os.listdir(ukbb_in_sample_ld_dir):
		if file_name.startswith('ukb_imp_v3_chimp.c' + chrom_num + '_') == False:
			continue
		if file_name.endswith('.compute_ld.sbatch.log'):
			continue
		full_file_name = ukbb_in_sample_ld_dir + file_name

		file_info = file_name.split('_')


		file_index_start = int(file_info[4].split('s')[1])
		file_index_end = int(file_info[5].split('e')[1])

		if indices_dont_lie_in_file(sample_ld_variant_indices, file_index_start, file_index_end):
			continue

		file_indices = []
		col_names = []
		for ii, sample_index in enumerate(sample_ld_variant_indices):
			if sample_index >= file_index_start and sample_index < file_index_end:
				file_indices.append(ii)
				col_names.append(sample_index - file_index_start)
		file_indices = np.asarray(file_indices)
		col_names = np.asarray(col_names)

		sparse_mat, sparse_mat_info = read_ld(full_file_name)
		ld_mat[:, file_indices] = (sparse_mat[sample_ld_variant_indices,:][:,col_names]).toarray()

	if np.sum(ld_mat == -2000.0) > 0:
		print('assumption eroror')
		pdb.set_trace()
	return ld_mat

def correct_ld_mat_for_af_standardization(ld_mat):
	n_snps = ld_mat.shape[0]
	correction = 1.0/np.diag(ld_mat)
	for snp_iter in range(n_snps):
		ld_mat[:,snp_iter] = ld_mat[:,snp_iter]*np.sqrt(correction[snp_iter])
		ld_mat[snp_iter,:] = ld_mat[snp_iter,:]*np.sqrt(correction[snp_iter])
	return ld_mat

def make_window_variant_info_file(window_variant_info_file, final_rsids, final_snp_pos, final_snp_alleles, chrom_num):
	t2 = open(window_variant_info_file,'w')
	t2.write('rsid\tchrom_num\tsnp_pos\tallele0\tallele1\n')

	for ii, rsid in enumerate(final_rsids):
		t2.write(rsid + '\t' + chrom_num + '\t' + str(final_snp_pos[ii]) + '\t' + final_snp_alleles[ii,0] + '\t' + final_snp_alleles[ii,1] + '\n')
	t2.close()
	return



chrom_num = sys.argv[1]
trait_summary_file = sys.argv[2]
overlapping_windows_file = sys.argv[3]
ukbb_in_sample_ld_dir = sys.argv[4]
ukbb_in_sample_genotype_dir = sys.argv[5]
processed_ukbb_in_sample_ld_dir = sys.argv[6]  # Output dir


# create mapping from RS_ID to UKBB in_sample LD variant INDEX
chrom_pvar_file = ukbb_in_sample_genotype_dir + 'ukb_imp_chr' + chrom_num + '_v3_chimp.pvar'
in_sample_rs_id_to_in_sample_variant, in_sample_rs_id_to_alleles, in_sample_rs_id_to_position, in_sample_rsids, in_sample_positions = create_mapping_from_rsid_to_in_sample_variant_index(chrom_pvar_file)


# Create dictionary list of summary statistics rs ids
sumstat_rsid_dictionary = create_summary_statistic_rsid_dictionary(trait_summary_file, chrom_num)

# Open output file
output_file = processed_ukbb_in_sample_ld_dir + 'in_sample_ld_summary_chrom_' + chrom_num + '.txt'
t = open(output_file,'w')

# Stream genome wide windows file
f = open(overlapping_windows_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\tvariant_info_file\tukbb_in_sample_ld_file\n')
		continue
	line_chrom_num = data[0]
	# Limit to windows on this chromosome
	if line_chrom_num != chrom_num:
		continue
	timer1 =time.time()
	# Info about this window
	window_name = data[0] + ':' + data[1] + ':' + data[2]
	print(window_name)
	# Window start and end
	window_start = int(data[1])
	window_end = int(data[2])

	# Extract in sample rsids
	window_indices = (in_sample_positions >= window_start) & (in_sample_positions < window_end)
	window_rsids = in_sample_rsids[window_indices]


	final_rsids = []
	final_snp_alleles = []
	final_snp_pos = []
	final_in_sample_flips = []
	final_in_sample_ld_indices = []

	for window_rsid in window_rsids:
		if window_rsid not in sumstat_rsid_dictionary:
			continue
		(sumstat_pos, sumstat_a0, sumstat_a1) = sumstat_rsid_dictionary[window_rsid]
		in_sample_pos = in_sample_rs_id_to_position[window_rsid]
		in_sample_a0, in_sample_a1 = in_sample_rs_id_to_alleles[window_rsid]

		# Error checking
		if sumstat_pos != in_sample_pos:
			print('snp position assumption error')
			pdb.set_trace()
		if sumstat_a0 != in_sample_a0 or sumstat_a1 != in_sample_a1:
			if sumstat_a0 != in_sample_a1 or sumstat_a1 != in_sample_a0:
				print('snp_allele assumption eroror')
				pdb.set_trace()

		# Get sign flips
		if sumstat_a0 == in_sample_a0 and sumstat_a1 == in_sample_a1:
			flipper = 1.0
		elif sumstat_a0 == in_sample_a1 and sumstat_a1 == in_sample_a0:
			flipper = -1.0
		else:
			print('assumption error')
			pdb.set_trace()

		final_rsids.append(window_rsid)
		final_snp_alleles.append((sumstat_a0, sumstat_a1))
		final_snp_pos.append(int(sumstat_pos))
		final_in_sample_flips.append(flipper)
		final_in_sample_ld_indices.append(in_sample_rs_id_to_in_sample_variant[window_rsid])

	# Convert to numpy arrays
	final_rsids = np.asarray(final_rsids)
	final_snp_alleles = np.asarray(final_snp_alleles)
	final_snp_pos = np.asarray(final_snp_pos)
	final_in_sample_flips = np.asarray(final_in_sample_flips)
	final_in_sample_ld_indices = np.asarray(final_in_sample_ld_indices)

	if len(final_rsids) <= 1:
		print(window_name + ' skipped do to having ' + str(len(final_rsids)) + ' snps in window')
		continue
	if np.array_equal(np.argsort(final_snp_pos), np.arange(len(final_snp_pos))) == False:
		print('positino sorting assumption eroror')
		pdb.set_trace()

	# NOW GET LD matrix
	ukbb_in_sample_ld_mat = extract_ld_mat_from_in_sample_ld(final_in_sample_ld_indices, ukbb_in_sample_ld_dir, chrom_num)
	# Flip alleles to make sure it matches summary statistics
	for var_index, sample_ld_flip_value in enumerate(final_in_sample_flips):
		if sample_ld_flip_value == -1.0:
			ukbb_in_sample_ld_mat[var_index, :] = ukbb_in_sample_ld_mat[var_index, :]*-1.0
			ukbb_in_sample_ld_mat[:, var_index] = ukbb_in_sample_ld_mat[:, var_index]*-1.0

	# Correct ld mat to make it a true correlation
	corrected_ukbb_in_sample_ld_mat = correct_ld_mat_for_af_standardization(ukbb_in_sample_ld_mat)

	# window variant info file
	window_variant_info_file = processed_ukbb_in_sample_ld_dir + window_name + '_variant_info.txt'
	make_window_variant_info_file(window_variant_info_file, final_rsids, final_snp_pos, final_snp_alleles, chrom_num)

	# window LD File
	window_ld_file = processed_ukbb_in_sample_ld_dir + window_name + '_in_sample_LD.npy'
	np.save(window_ld_file, corrected_ukbb_in_sample_ld_mat)

	# Print to output
	t.write(line + '\t' + window_variant_info_file + '\t' + window_ld_file + '\n')
	timer2=time.time()
	print(timer2-timer1)
	
f.close()
t.close()
