# %%
import re
import os
import sys
import glob
import math
import pysam
import argparse
import pybedtools
import subprocess
import numpy as np
import polars as pl
from pyfaidx import Fasta
from functools import partial
from datetime import datetime

def filter_reads(segment):
	if segment.mapping_quality > 0:
		return True
	else:
		return False

def output_reads_alin_info(sv_id, ref_bam, pse_bam, alt_sv_info, chrom, sv_start, sv_end, ref_chr_len, alt_chr_len, slop_len = 150):
	read_ids = []
	read_lens = []
	sv_ids = []
	allel_ids = []
	cigars  = []
	NMs = []
	match_bases = []
	ref_starts = []
	reads_from = []

	start = max(sv_start - slop_len, 0) 
	end = min(sv_end + slop_len, ref_chr_len)

	clip_line = 3

	for segment in ref_bam.fetch(chrom, start, end):
		if segment.mapping_quality == 0 or segment.cigartuples is None:
			continue

		query_name = segment.query_name
		query_seq = segment.query_sequence
		query_len = len(query_seq)

		s = segment.reference_start
		e = segment.reference_end
		match_base_count = sum([tup[1] for tup in segment.cigartuples if tup[0] == 0])
		nm = segment.get_tag("NM")
		
		readsuffix = '_1' if segment.is_read1 else '_2'

		alt_support = False
		if segment.cigarstring and re.search("[SH]", segment.cigarstring):
			no_num_cigarstr = re.sub(r'\d+', '', segment.cigarstring)
			if re.search(r'^[SH].*[SH]$', no_num_cigarstr):
				if -clip_line < s - sv_end < clip_line:
					alt_support = True
				elif -clip_line < e - sv_start < clip_line:
					alt_support = True
			elif (re.search(r'^[SH]', no_num_cigarstr) and -clip_line < s - sv_end < clip_line):
				alt_support = True
			elif (re.search(r'[SH]$', no_num_cigarstr) and -clip_line < e - sv_start < clip_line):
				alt_support = True

		elif nm < 3 and match_base_count == query_len and not segment.is_secondary and not re.search("[SH]", segment.cigarstring):
			if ((sv_start - query_len) < s < sv_start or (sv_end - query_len) < s < sv_end):
					read_ids.append(f'{query_name}{readsuffix}')
					read_lens.append(query_len)
					sv_ids.append(sv_id)
					allel_ids.append(0)
					cigars.append(segment.cigarstring)
					NMs.append(nm)
					match_bases.append(match_base_count)
					ref_starts.append(s)
					reads_from.append(0)

		if alt_support == True:
			read_ids.append(f'{query_name}{readsuffix}_clip')
			read_lens.append(query_len)
			sv_ids.append(sv_id)
			allel_ids.append(1)
			cigars.append(segment.cigarstring)
			NMs.append(0)
			match_bases.append(slop_len)
			ref_starts.append(segment.reference_start)
			reads_from.append(0)

	pse_chr, pse_sv_start, pse_sv_end = alt_sv_info[sv_id]
	pse_start = max(pse_sv_start - slop_len, 0)
	pse_end = min(pse_sv_end + slop_len, alt_chr_len)

	for segment in pse_bam.fetch(pse_chr, pse_start, pse_end):
		if segment.mapping_quality == 0 or segment.cigartuples is None:
			continue
		
		query_name = segment.query_name
		query_seq = segment.query_sequence
		query_len = len(query_seq)

		s = segment.reference_start
		e = segment.reference_end
		match_base_count = sum([tup[1] for tup in segment.cigartuples if tup[0] == 0])
		nm = segment.get_tag("NM")

		readsuffix = '_1' if segment.is_read1 else '_2'

		alt_support = False
		if segment.cigarstring and re.search("[SH]", segment.cigarstring):
			no_num_cigarstr = re.sub(r'\d+', '', segment.cigarstring)
			if re.search(r'^[SH].*[SH]$', no_num_cigarstr):
				if -clip_line < s - pse_sv_end < clip_line:
					alt_support = True
				elif -clip_line < e - pse_sv_start < clip_line:
					alt_support = True
			elif (re.search(r'^[SH]', no_num_cigarstr) and -clip_line < s - pse_sv_end < clip_line):
				alt_support = True
			elif (re.search(r'[SH]$', no_num_cigarstr) and -clip_line < e - pse_sv_start < clip_line):
				alt_support = True
		elif nm < 3 and match_base_count == query_len and not segment.is_secondary and not re.search("[SH]", segment.cigarstring):
			if ((pse_sv_start - query_len) < s < pse_sv_start or (pse_sv_end - query_len) < s < pse_sv_end):
					read_ids.append(f'{query_name}{readsuffix}')
					read_lens.append(query_len)
					sv_ids.append(sv_id)
					allel_ids.append(1)
					cigars.append(segment.cigarstring)
					NMs.append(nm)
					match_bases.append(match_base_count)
					ref_starts.append(s)
					reads_from.append(1)

		if alt_support == True:
			read_ids.append(f'{query_name}{readsuffix}_clip')
			read_lens.append(query_len)
			sv_ids.append(sv_id)
			allel_ids.append(0)
			cigars.append(segment.cigarstring)
			NMs.append(0)
			match_bases.append(slop_len)
			ref_starts.append(segment.reference_start)
			reads_from.append(1)


	df = pl.DataFrame({'read_id': read_ids,
					   'read_len': read_lens,
					   'sv_id': sv_ids,
					   'allel_id': allel_ids,
					   'r_start': ref_starts,
					   'cigar': cigars,
					   'NM': NMs,
					   'match_base': match_bases,
					   'read_from': reads_from})
	
	return df

def reads_best_scoring(read_info_df):
	df = (read_info_df.lazy()
		  .with_columns(
			  (pl.col("match_base") - (pl.col('NM') * 2)).alias('read_score'))
		  .with_columns(
			  pl.col("read_score").max().over("read_id").alias("max_score_per_read"))
		  .filter(pl.col("read_score") == pl.col("max_score_per_read"))
		  .with_columns(
			  pl.col("read_id").count().over("read_id").alias("num_besthits_per_read"))
		  .unique()
		 ).collect()
	
	df = df.sort(['read_id', 'sv_id', 'allel_id'])

	return df

def get_read_counts(best_reads_file):
	gdf = (best_reads_file.lazy()
		   .filter(pl.col('num_besthits_per_read') == 1)
		   .group_by(["sv_id", "allel_id"])
		   .agg([pl.count("read_id").alias('read_counts')])
		   .sort(['sv_id', 'allel_id', 'read_counts'], descending=[False, False, True])
		  ).collect()
	
	return gdf

def log_choose(n, k):
	r = 0.0

	if k * 2 > n:
		k = n - k

	for d in range(1,k+1):
		r += math.log(n, 10)
		r -= math.log(d, 10)
		n -= 1

	return r

def bayes_gt(ref, alt, is_depth=False):
	if is_depth:
		p_alt = [0.0625, 0.5, 0.99]
	else:
		p_alt = [1e-3, 0.5, 0.9]

	total = ref + alt
	log_combo = log_choose(total, alt)

	lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
	lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
	lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)
	
	gt_types = ('0/0', '0/1', '1/1')
	pre_gt = gt_types[np.argmax((lp_homref, lp_het, lp_homalt))]
	
	return [pre_gt, lp_homref, lp_het, lp_homalt]

def bayes_gt_wrapper(index, x, y, is_depth=False):
	values = bayes_gt(int(x), int(y), is_depth)
	return values[index]

bayes_gt_1 = partial(bayes_gt_wrapper, 0)
bayes_gt_2 = partial(bayes_gt_wrapper, 1)
bayes_gt_3 = partial(bayes_gt_wrapper, 2)
bayes_gt_4 = partial(bayes_gt_wrapper, 3)

def get_all_ref_depth(script_dir, ref_bam_file, out_file_prefix, ref_fasta_file, thread=1):
	software_path = os.path.join(script_dir, 'mosdepth')
	
	mosdepth_command = [
		software_path,
		'--threads', str(thread),
		'--fast-mode', '--no-per-base',
		'--by', "300",
		f'{out_file_prefix}.mosdepth', ref_bam_file]
	
	try:
		subprocess.run(mosdepth_command, check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error occurred while running mosdepth: {e}")
		sys.exit(1)
	except Exception as e:
		print(f"Error occurred while running mosdepth: {e}")
		sys.exit(1)
		
	out_mosdepth_file = f'{out_file_prefix}.mosdepth.regions.bed.gz'
	ref_fasta = Fasta(ref_fasta_file)
	ref_depth_bed = pybedtools.BedTool(out_mosdepth_file)
	depth_gc_lines = []
	for window in ref_depth_bed:
		gc = ref_fasta[window.chrom][window.start:window.end].gc
		depth_gc_lines.append([window.chrom, window.start, window.end, float(window.fields[3]), gc])

	columns_name = ['chr', 'start', 'end', 'depth', 'gc']
	depth_gc_df = pl.DataFrame(depth_gc_lines) 
	depth_gc_df = pl.DataFrame({new_col: depth_gc_df[old_col] for old_col, new_col in zip(depth_gc_df.columns, columns_name)}
							   ).with_columns((pl.col('gc') * 100).round().cast(pl.Int64).alias('gc(%)'))
	
	depth_GC_dict = {}
	for i in range(0, 101, 1):
		df_dep = depth_gc_df.filter(pl.col('gc(%)') == i)['depth']
		RD_median = df_dep.median()
		depth_GC_dict[i] = RD_median

	corr_40 = depth_GC_dict[40]
	depth_gc_df = (depth_gc_df.lazy()
				   .with_columns(pl.col('gc(%)')
								 .map_elements(lambda gc_content: depth_GC_dict.get(gc_content, None)).alias('median_depth'))
				   .with_columns((pl.col('depth') * (corr_40 / pl.col('median_depth'))).round(3).alias('corrected_depth'))
				   .fill_nan(0)
				   .fill_null(0)
				   ).collect()
	
	all_ref_median_depth = depth_gc_df['corrected_depth'].median()
	
	return depth_GC_dict, all_ref_median_depth


def get_output_file_path(output_dir_path, file_name):
	"""
	构建输出文件的完整路径
	Args:
		output_dir_path: 输出目录的路径
		file_name: 输出文件名
	Returns:
		输出文件的完整路径
	"""
	return os.path.join(output_dir_path, file_name)

def main(args=None):
	parser = argparse.ArgumentParser(description="This script gets breakpoint reads signatures for SV sites from ref and alt bam.")

	parser.add_argument("--ref_fasta", type=str, required=True, help="The ref genome Fasta file", metavar="file")
	parser.add_argument("--alt_fasta", type=str, required=True, help="The alt genome Fasta file", metavar="file")
	parser.add_argument("--ref_sv_vcf", type=str, required=True, help="The VCF file based ref genome", metavar="file")
	parser.add_argument("--alt_sv_bed", type=str, required=True, help="The BED file based alt genome", metavar="file")
	parser.add_argument("--ref_bam", type=str, required=True, help="The ref genome bam file", metavar="file")
	parser.add_argument("--alt_bam", type=str, required=True, help="The alt genome bam file", metavar="file")
	parser.add_argument("--read_leagth", type=str, default="150", help="The read length. Default is 150. optional:{100, 150, 250}", metavar="int")
	parser.add_argument("-t", "--threads", type=str, default="1", help="threads. No more than 4", metavar="int")
	parser.add_argument("-o", "--out", type=str, default="BreakPoint_ReadDepth_feature_out", help="The output dir name. Default: BreakPoint_ReadDepth_feature_out", metavar="dir")

	parsed_args = parser.parse_args(args=args)
	ref_vcf_file = parsed_args.ref_sv_vcf
	alt_bed_file = parsed_args.alt_sv_bed
	ref_bam_file = parsed_args.ref_bam
	alt_bam_file = parsed_args.alt_bam
	ref_fasta_file = parsed_args.ref_fasta
	alt_fasta_file = parsed_args.alt_fasta
	read_length = parsed_args.read_leagth
	threads = parsed_args.threads
	output_dir = parsed_args.out

	threads = int(threads)
	slop_length = int(read_length)
	output_dir_path = os.path.abspath(output_dir)
	os.makedirs(output_dir_path, exist_ok=True)

	alt_sv_info = {}
	with open(alt_bed_file, "r") as f:
		for record in f:
			items = record.strip().split('\t')
			sv_id = items[3]
			chr = str(items[0])
			pos = int(items[1])
			end = int(items[2])
			alt_sv_info[sv_id] = [chr, pos, end]

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "Read vcf file Finished!", "There are " + str(len(alt_sv_info)) + " SVs.")

	script_dir = os.path.dirname(os.path.abspath(__file__))
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "mosdepth: ", script_dir+"/mosdepth")
	out_file_prefix = get_output_file_path(output_dir_path, "tmp")
	depth_GC_dict, ref_median_depth = get_all_ref_depth(script_dir, ref_bam_file, out_file_prefix, ref_fasta_file, threads)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "Global median depth in ref bam: ", str(ref_median_depth))

	bp_out_file = get_output_file_path(output_dir_path, "BP_2Bam_result.tsv")

	with open(bp_out_file, 'w') as output_file:
		output_file.write("sv_id\t0\t1\tread_counts_sum\tBP_GT_pred\tBP_lp_homref\tBP_lp_het\tBP_lp_homalt\n")

	all_read_counts_df = pl.DataFrame()
	ref_bam = pysam.AlignmentFile(ref_bam_file, "r")
	pse_bam = pysam.AlignmentFile(alt_bam_file, "r")
	ref_sv_vcf = pysam.VariantFile(ref_vcf_file)
	ref_genome = Fasta(ref_fasta_file)
	alt_genome = Fasta(alt_fasta_file)

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "Extracting features ...")
	read_depth_out_list = []
	correct_40 = depth_GC_dict[40]
	with open(bp_out_file, 'a') as output_file:
		for record in ref_sv_vcf:
			chrom = record.chrom
			sv_id = record.id
			ref_seq = record.ref
			ref_start = record.start
			ref_end = len(ref_seq) + ref_start
			sv_type = record.info['SVTYPE']
			ref_chr_len = len(ref_genome[chrom])
			alt_chr_len = len(alt_genome[chrom])
			alt_chrom, alt_start, alt_end = alt_sv_info[sv_id]

			read_info_df = output_reads_alin_info(
				sv_id, ref_bam, pse_bam, alt_sv_info, chrom, ref_start, ref_end, ref_chr_len, alt_chr_len, slop_length)
			
			reads_best_df = reads_best_scoring(read_info_df)
 
			read_counts_df = get_read_counts(reads_best_df)

			if len(read_counts_df) == 0:
				outitems = [sv_id, '-', '-', '-', './.','-','-','-']
				outline = '\t'.join([str(x) for x in outitems]) + '\n'
				output_file.write(outline)
			else:
				all_read_counts_df = pl.concat([all_read_counts_df, read_counts_df])

			if sv_type == "DEL":
				gc = ref_genome[chrom][ref_start:ref_end].gc
				depths_region = ref_bam.count_coverage(chrom, ref_start, ref_end, quality_threshold = 20, read_callback=filter_reads)
				depth_list = np.sum(depths_region, axis=0)
				depth = np.median(depth_list)
			elif sv_type == "INS":
				gc = alt_genome[chrom][alt_start:alt_end].gc
				depths_region = pse_bam.count_coverage(chrom, alt_start, alt_end, quality_threshold = 20, read_callback=filter_reads)
				depth_list = np.sum(depths_region, axis=0)
				depth = np.median(depth_list)
			
			ref_front_1k = max(ref_start - 1000, 0)
			ref_end_1k = min(ref_end + 1000, ref_chr_len)
			alt_front_1k = max(alt_start - 1000, 0)
			alt_end_1k = min(alt_end + 1000, alt_chr_len)

			ref_front_1k_depth = np.median(np.sum(ref_bam.count_coverage(
				chrom, ref_front_1k, ref_start, quality_threshold = 20, read_callback=filter_reads), axis=0))
			ref_end_1k_depth = np.median(np.sum(ref_bam.count_coverage(
				chrom, ref_end, ref_end_1k, quality_threshold = 20, read_callback=filter_reads), axis=0))
			alt_front_1k_depth = np.median(np.sum(pse_bam.count_coverage(
				chrom, alt_front_1k, alt_start, quality_threshold = 20, read_callback=filter_reads), axis=0))
			alt_end_1k_depth = np.median(np.sum(pse_bam.count_coverage(
				chrom, alt_end, alt_end_1k, quality_threshold = 20, read_callback=filter_reads), axis=0))

			if (ref_front_1k_depth + alt_front_1k_depth) > 0:
				front_depth = ref_front_1k_depth / (ref_front_1k_depth + alt_front_1k_depth)
			else:
				front_depth = 0.5

			if (ref_end_1k_depth + alt_end_1k_depth) > 0:
				end_depth = ref_end_1k_depth / (ref_end_1k_depth + alt_end_1k_depth)
			else:
				end_depth = 0.5

			gc_percentage = int(gc * 100)

			if depth_GC_dict[gc_percentage] != 0 and depth_GC_dict[gc_percentage] is not None:
				depth_corr = round((correct_40 / float(depth_GC_dict[gc_percentage])) * depth, 3)
			else:
				depth_corr = depth
				
			read_depth_out_list.append({
				'sv_id': sv_id,
				'sv_type': sv_type,
				'sv_RD': depth,
				'corrected_sv_RD': depth_corr,
				'sv_gc_content': round(gc, 3),
				'front_depth': front_depth,
				'end_depth': end_depth
				})

	bp_genotype_df = (
		all_read_counts_df.pivot(index="sv_id", columns="allel_id",values="read_counts", aggregate_function=None).fill_null(0)
		.with_columns((pl.col("0") + pl.col("1")).alias("BP_read_counts_sum"))
		.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_1(x['0'], x['1'])).alias('BP_GT_pred')))
		.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_2(x['0'], x['1'])).alias('BP_lp_homref')))
		.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_3(x['0'], x['1'])).alias('BP_lp_het')))
		.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_4(x['0'], x['1'])).alias('BP_lp_homalt')))
		)

	rd_genotype_df = pl.DataFrame(read_depth_out_list)
	rd_genotype_df = rd_genotype_df.with_columns(
		pl.when(pl.col('sv_type') == 'DEL')
		.then(pl.col('corrected_sv_RD'))
		.when(pl.col('sv_type') == 'INS').then(
			pl.when(pl.col('corrected_sv_RD') <= ref_median_depth)
			.then(ref_median_depth - pl.col('corrected_sv_RD'))
			.otherwise(0))
		.otherwise(-1)
		.alias("0"))

	rd_genotype_df = rd_genotype_df.with_columns(
		pl.when(pl.col('sv_type') == 'INS')
		.then(pl.col('corrected_sv_RD'))
		.when(pl.col('sv_type') == 'DEL')
		.then(
			pl.when(pl.col('corrected_sv_RD') <= ref_median_depth)
			.then(ref_median_depth - pl.col('corrected_sv_RD'))
			.otherwise(0))
		.otherwise(-1)
		.alias("1"))
		
	rd_genotype_df = (rd_genotype_df
					.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_1(x['0'], x['1'], is_depth=True)).alias('RD_GT_pred')))
					.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_2(x['0'], x['1'], is_depth=True)).alias('RD_lp_homref')))
					.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_3(x['0'], x['1'], is_depth=True)).alias('RD_lp_het')))
					.with_columns((pl.struct('0', '1').map_elements(lambda x: bayes_gt_4(x['0'], x['1'], is_depth=True)).alias('RD_lp_homalt')))
					)

	rd_out_file = get_output_file_path(output_dir_path, "RD_2Bam_result.tsv")
	rd_genotype_df.write_csv(rd_out_file, include_header=True, separator="\t")


	with open(bp_out_file, 'a') as output_file:
		for row in bp_genotype_df.rows():
			row_str = "\t".join(map(str, row)) + "\n"
			output_file.write(row_str)


	ref_bam.close()
	pse_bam.close()
	ref_sv_vcf.close()
	ref_genome.close()
	alt_genome.close()

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "BreakPoint and ReadDepth feature extraction is complete!")

	df_bp = pl.read_csv(bp_out_file,  separator="\t").select(["sv_id", "BP_lp_homref", "BP_lp_het", "BP_lp_homalt"])
	df_rd = rd_genotype_df.select(["sv_id", "sv_gc_content", "front_depth","end_depth", "RD_lp_homref", "RD_lp_het", "RD_lp_homalt"])

	sv_feature_df = df_rd.join(df_bp, on="sv_id", how="inner")
	featrue_out = get_output_file_path(output_dir_path, "BreakPoint_ReadDepth_2Bam_feature.tsv")
	sv_feature_df.write_csv(featrue_out, include_header=True, separator="\t")
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "BreakPoint and ReadDepth features are extracted and saved to file:", featrue_out)

	featrue_out = get_output_file_path(output_dir_path, "*.mosdepth*")
	mosdepth_files = glob.glob(featrue_out)
	[os.remove(file_path) for file_path in mosdepth_files]
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "DONE!")


if __name__ == "__main__":
	main()
