import argparse
import pysam
import numpy as np
import polars as pl
from datetime import datetime

def current_time():
	return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def main(args=None):
	parser = argparse.ArgumentParser(description="integrate paragraph results para_feature.tsv")
	parser.add_argument("--ref_paragraph_vcf", type=str, required=True, help="paragraph genotyped VCF based ref genome", metavar="file")
	parser.add_argument("--alt_paragraph_vcf", type=str, required=True, help="paragraph genotyped VCF based alt genome", metavar="file")
	parser.add_argument("--ref_sample_id", type=str, default="GT_paragraph_ref", help="sample ID in ref paragraph VCF", metavar="file")
	parser.add_argument("--alt_sample_id", type=str, default="GT_paragraph_alt", help="sample ID in alt paragraph VCF", metavar="file")
	parser.add_argument("-o", "--out", type=str, default="para_feature.tsv", help="output paragraph feature matrix", metavar='file')

	parsed_args = parser.parse_args(args=args)
	ref_paragraph_vcf_file = parsed_args.ref_paragraph_vcf
	alt_paragraph_vcf_file = parsed_args.alt_paragraph_vcf
	ref_sample_id = parsed_args.ref_sample_id
	alt_sample_id = parsed_args.alt_sample_id
	output_file = parsed_args.out

	ref_gt_dict = {
		(0, 1): 1,
		(1, 0): 1,
		(0, 0): 0,
		(1, 1): 2,
		(0,): 0,
		(1,): 2,
		(None,): np.nan,
	}

	alt_gt_dict = {
		(0, 1): 1,
		(1, 0): 1,
		(0, 0): 2,
		(1, 1): 0,
		(0,): 2,
		(1,): 0,
		(None,): np.nan,
	}

	para_FT = {
		"PASS": 0,
		"BP_DEPTH": 1,
		"BP_NO_GT": 2, 
		"CONFLICT": 3,
		"GQ": 4,
		('UNMATCHED', 'NO_VALID_GT'): 5,
		('NO_VALID_GT', 'UNMATCHED'): 5,
		('GQ', 'BP_DEPTH'): 6,
		('BP_DEPTH', 'GQ'): 6,
	}

	ref_paragraph_gt_info = []
	ref_paragraph_vcf = pysam.VariantFile(ref_paragraph_vcf_file)
	for record in ref_paragraph_vcf:
		sv_id = record.id
		ref_para_geno_info = record.samples.get(ref_sample_id, {})
		ref_sv_gt_tup = ref_para_geno_info.get('GT', (None,))
		if ref_sv_gt_tup == (None,):
			ref_sv_PL = np.nan
		else:
			ref_sv_PL = float(sorted(ref_para_geno_info['PL'])[1])
		
		ref_sv_FT_tup = ref_para_geno_info['FT']
		ref_sv_DP = ref_para_geno_info['DP']

		ref_sv_gt = float(ref_gt_dict.get(ref_sv_gt_tup, np.nan))
		ref_sv_FT = para_FT.get(ref_sv_FT_tup, np.nan)

		ref_paragraph_gt_info.append({
			'sv_id': sv_id,
			'ref_GT': ref_sv_gt,
			'ref_FT': ref_sv_FT,
			'ref_DP': ref_sv_DP,
			'ref_PL': ref_sv_PL
		})
	
	ref_info_df = pl.DataFrame(ref_paragraph_gt_info)
	ref_paragraph_vcf.close()


	alt_paragraph_gt_info = []
	alt_paragraph_vcf = pysam.VariantFile(alt_paragraph_vcf_file)
	for record in alt_paragraph_vcf:
		sv_id = record.id
		alt_para_geno_info = record.samples.get(alt_sample_id, {})
		alt_sv_gt_tup = alt_para_geno_info.get('GT', (None,))
		if alt_sv_gt_tup == (None,):
			alt_sv_PL = np.nan
		else:
			alt_sv_PL = float(sorted(alt_para_geno_info['PL'])[1])

		alt_sv_FT_tup = alt_para_geno_info['FT']
		alt_sv_DP = alt_para_geno_info['DP']

		alt_sv_gt = float(alt_gt_dict.get(alt_sv_gt_tup, np.nan))
		alt_sv_FT = para_FT.get(alt_sv_FT_tup, np.nan)

		alt_paragraph_gt_info.append({
			'sv_id': sv_id,
			'alt_GT': alt_sv_gt,
			'alt_FT': alt_sv_FT,
			'alt_DP': alt_sv_DP,
			'alt_PL': alt_sv_PL
		})

	alt_info_df = pl.DataFrame(alt_paragraph_gt_info)

	para_info_df = ref_info_df.join(alt_info_df, on="sv_id", how="left")

	para_info_df = (
		para_info_df
		.with_columns((pl.col("ref_DP") - pl.col("alt_DP")).alias("DP"))
		.with_columns((pl.col("ref_PL") - pl.col("alt_PL")).alias("PL"))
	)

	selected_columns = ['sv_id', 'ref_GT', 'alt_GT', 'ref_FT', 'alt_FT', 'DP', 'PL']
	para_info_df.select(selected_columns).write_csv(output_file, separator="\t")

	print(f"{current_time()} Successful generation {output_file}")
	print("ALL DONE.")

if __name__ == "__main__":
	main()
