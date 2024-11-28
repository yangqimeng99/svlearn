# %%
import os
import pysam
import argparse
import pybedtools
import numpy as np
import pandas as pd
from functools import reduce
from datetime import datetime

def process_TRF_out(trf_file):
	trf_df = pd.read_csv(trf_file, sep='\t|;|=', engine='python',
					header=None, 
					usecols=[0,3,4,11,13], 
					names=['chr', 'start', 'end', 'period', 'copies'],
					dtype={'chr': str, 'start': int, 'end': int, 'period': int, 'copies': float})

	trf_df['start'] = trf_df['start'] - 1
	trf_df['type'] = trf_df['period'].apply(lambda x: 'STR' if x <= 6 else 'VNTR')
	trf_VNTR_bed = pybedtools.BedTool.from_dataframe(trf_df[trf_df['type'] == 'VNTR']).sort().merge()
	trf_STR_bed = pybedtools.BedTool.from_dataframe(trf_df[trf_df['type'] == 'STR']).sort().merge()
	trf_bed = pybedtools.BedTool.from_dataframe(trf_df).sort().merge()

	return trf_bed, trf_VNTR_bed, trf_STR_bed


def process_vcf_to_del_bed(ref_vcf):
	variants = []
	for record in ref_vcf:
		sv_type = record.info.get("SVTYPE")
		if sv_type == "DEL":
			sv_id = record.id
			chrom = record.chrom
			start = record.start
			end = start + len(record.ref)
			interval = f"{chrom}\t{start}\t{end}\t{sv_id}\t{sv_type}"
			variants.append(interval)

	df_dels = pd.DataFrame([x.split('\t') for x in variants], columns=['chrom', 'start', 'end', 'sv_id', 'sv_type'])
	bed_dels = pybedtools.BedTool.from_dataframe(df_dels)

	return bed_dels


def filter_del(feature):
	return "DEL" not in feature.fields[4]


def process_repeat_masker(data, repeat_type):
	filtered_data = data[data[10].str.contains(repeat_type)].copy()
	filtered_data.loc[:, 5] = filtered_data[5].astype(int)
	filtered_data.loc[:, 6] = filtered_data[6].astype(int)
	filtered_data.loc[:, 5] = filtered_data[5] - 1
	
	filtered_data_bed = filtered_data[[4, 5, 6]]
	
	bed = pybedtools.BedTool.from_dataframe(filtered_data_bed)
	bed_sorted = bed.sort()
	bed_merged = bed_sorted.merge()

	return bed_merged
	

# %%
def process_biser_out(biser_file, Satellite_bed):
	dtype_spec = {
	0: str, 1: int, 2: int, 3: str, 4: int, 5: int, 
	6: str, 7: float, 8: str, 9: str, 10: int, 11: int, 12: str, 13: str
	}

	df = pd.read_csv(biser_file, sep='\t', header=None, dtype=dtype_spec)
	df.drop(columns=[6, 12], inplace=True)

	df[['x', 'x_num', 'ID', 'ID_num']] = df.iloc[:, -1].str.split('[;=]', expand=True)
	df[['x_num', 'ID_num']] = df[['x_num', 'ID_num']].astype(float)
	
	df = df[(df['ID_num'] <= 50) & (df[11] >= 1000) & (df['x_num'] < 10)] 
	df.drop(df.columns[-5], axis=1, inplace=True)

	bed_a = pybedtools.BedTool.from_dataframe(df[[0, 1, 2, 3, 4, 5]])

	intersect_result = bed_a.intersect(Satellite_bed, wo=True).sort().groupby(g=[1,2,3,4,5,6], c=10, o=['sum']).to_dataframe()
	intersect_result['overlap_ratio'] = intersect_result['thickStart'] / (intersect_result['end'] - intersect_result['start'])
	filtered_overlap = intersect_result[intersect_result['overlap_ratio'] > 0.7]
	filtered_overlap = filtered_overlap.astype({'chrom': 'str', 'name': 'str'})
	merged_df = pd.merge(df, filtered_overlap, how='left', indicator=True,
						left_on=[0, 1, 2, 3, 4, 5],
						right_on=['chrom', 'start', 'end', 'name', 'score', 'strand'])
	merged_df_filted = merged_df[merged_df['_merge'] == 'left_only']

	merged_df_filted_1 = merged_df_filted.iloc[:, :3]
	merged_df_filted_2 = merged_df_filted.iloc[:, 3:6]
	merged_df_filted_2.columns = merged_df_filted_1.columns

	df_combined = pd.concat([merged_df_filted_1, merged_df_filted_2], ignore_index=True)
	df_sorted = df_combined.sort_values(by=[0, 1]).reset_index(drop=True)
	df_combined_bed_merge = pybedtools.BedTool.from_dataframe(df_sorted).merge()

	return df_combined_bed_merge


def sv_class_based_SDs(sv_bed, sd_bed):
	svBed_intersect_SDsBed = (sv_bed.intersect(sd_bed, wao=True)
						   .to_dataframe(names=['chr', 'start', 'end', 'sv_id', 'sv_type', 'sd_chr', 'sd_start', 'sd_end', 'overlap_length'], dtype={'chr': str}))

	svBed_intersect_SDsBed_length = (svBed_intersect_SDsBed[['sv_id', 'overlap_length']]
								   .groupby('sv_id')['overlap_length'].sum().reset_index())
	
	svBed_df = sv_bed.to_dataframe(names=['chr', 'start', 'end', 'sv_id', 'sv_type'])

	df_merged = pd.merge(svBed_df, svBed_intersect_SDsBed_length, on='sv_id', how='left')
	df_merged['sv_length'] = df_merged['end'] - df_merged['start']
	df_merged['SDs_content'] = df_merged['overlap_length'] / df_merged['sv_length']

	df_merged['SDs_class'] = df_merged.apply(lambda x: 1 if (x['SDs_content'] >= 0.8 or x['overlap_length'] >= 200) else 0, axis=1)
	sv_class_sd = df_merged[['sv_id', 'SDs_content', 'SDs_class']]

	return sv_class_sd


def svBed_intersect_Repeat(sv_bed, repeat_bed, repeat_type):
	svBed_intersect_repeatsBed = (sv_bed.intersect(repeat_bed, wao=True)
							  .sort()
							  .groupby(g=[1,2,3,4,5], c=9, o=['sum'])
							  .to_dataframe(names=['chr', 'start', 'end', 'sv_id', 'sv_type' ,repeat_type], dtype={'chr': str}))
	
	svBed_intersect_repeatsBed['sv_length'] = svBed_intersect_repeatsBed['end'] - svBed_intersect_repeatsBed['start']
	
	columns = list(svBed_intersect_repeatsBed.columns)
	columns[-2], columns[-1] = columns[-1], columns[-2]
	svBed_intersect_repeatsBed = svBed_intersect_repeatsBed[columns]

	return svBed_intersect_repeatsBed


def merge_dfs(left, right):
	key = ['chr', 'start', 'end', 'sv_id', 'sv_type', 'sv_length']
	return pd.merge(left, right, on=key, how='outer')



def calculate_TE_class_content(row):
	if row.iloc[6] / row.iloc[5] >= 0.8:
		return pd.Series([1, row.iloc[6] / row.iloc[5]])
	elif row.iloc[7] / row.iloc[5] >= 0.8:
		return pd.Series([2, row.iloc[7] / row.iloc[5]])
	elif row.iloc[8] / row.iloc[5] >= 0.8:
		return pd.Series([3, row.iloc[8] / row.iloc[5]])
	elif row.iloc[9] / row.iloc[5] >= 0.8:
		return pd.Series([4, row.iloc[9] / row.iloc[5]])
	elif (row.iloc[7] + row.iloc[8] + row.iloc[9] + row.iloc[6]) / row.iloc[5] >= 0.8:
		return pd.Series([5, (row.iloc[7] + row.iloc[8] + row.iloc[9] + row.iloc[6]) / row.iloc[5]])
	else:
		return pd.Series([0, (row.iloc[7] + row.iloc[8] + row.iloc[9] + row.iloc[6]) / row.iloc[5]])

def calculate_Satellite_class_content(row):
	if row.iloc[6] / row.iloc[5] >= 0.8:
		return pd.Series([6, row.iloc[6] / row.iloc[5]])
	else:
		return pd.Series([0, row.iloc[6] / row.iloc[5]])

def calculate_TRF_class_content(row):
	if row.iloc[6] / row.iloc[5] >= 0.8:
		return pd.Series([7, row.iloc[6] / row.iloc[5]])
	elif row.iloc[7] / row.iloc[5] >= 0.8:
		return pd.Series([8, row.iloc[7] / row.iloc[5]])
	elif (row.iloc[6] + row.iloc[7]) / row.iloc[5] >= 0.8:
		return pd.Series([9,(row.iloc[6] + row.iloc[7]) / row.iloc[5]])
	else:
		return pd.Series([0, (row.iloc[6] + row.iloc[7]) / row.iloc[5]])
	

def sv_class_based_Repeats(sv_bed, RM_bed_dict, trf_VNTR_bed, trf_STR_bed):

	columns_out = ['chr', 'start', 'end', 'sv_id', 'sv_type', 'sv_length', 'repeat_class', 'repeat_content']

	# RepeatMasker TE
	sv_intersect_LINE_df = svBed_intersect_Repeat(sv_bed, RM_bed_dict['LINE'], 'LINE')
	sv_intersect_SINE_df = svBed_intersect_Repeat(sv_bed, RM_bed_dict['SINE'], 'SINE')
	sv_intersect_LTR_df = svBed_intersect_Repeat(sv_bed, RM_bed_dict['LTR'], 'LTR')
	sv_intersect_DNA_df = svBed_intersect_Repeat(sv_bed, RM_bed_dict['DNA'], 'DNA')
	list_merge_TE_df = [sv_intersect_LINE_df, sv_intersect_SINE_df, sv_intersect_LTR_df, sv_intersect_DNA_df]
	TE_class_content_df = reduce(merge_dfs, list_merge_TE_df)
	TE_class_content_df[['repeat_class', 'repeat_content']] = TE_class_content_df.apply(calculate_TE_class_content, axis=1)
	TE_class_content_df = TE_class_content_df[columns_out]
	TEonly_class_content_df = TE_class_content_df[TE_class_content_df['repeat_class'] != 0]

	# Satellite
	filtered_TE_df = TE_class_content_df[TE_class_content_df['repeat_class'] == 0]
	filtered_TE_df = filtered_TE_df.iloc[:, :5]
	filtered_TE_bed = pybedtools.BedTool.from_dataframe(filtered_TE_df)
	Satellite_class_content_df = svBed_intersect_Repeat(filtered_TE_bed, RM_bed_dict['Satellite'], 'Satellite')
	Satellite_class_content_df[['repeat_class', 'repeat_content']] = Satellite_class_content_df.apply(calculate_Satellite_class_content, axis=1)
	Satellite_class_content_df = Satellite_class_content_df[columns_out]
	Satelliteonly_class_content_df = Satellite_class_content_df[Satellite_class_content_df['repeat_class'] != 0]

	# TRF
	filtered_TE_Satellite_df = Satellite_class_content_df[Satellite_class_content_df['repeat_class'] == 0]
	filtered_TE_Satellite_df = filtered_TE_Satellite_df.iloc[:, :5]
	filtered_TE_Satellite_bed = pybedtools.BedTool.from_dataframe(filtered_TE_Satellite_df)
	sv_intersect_VNTR_df = svBed_intersect_Repeat(filtered_TE_Satellite_bed, trf_VNTR_bed, 'VNTR')
	sv_intersect_STR_df = svBed_intersect_Repeat(filtered_TE_Satellite_bed, trf_STR_bed, 'STR')
	list_merge_TRF_df = [sv_intersect_VNTR_df, sv_intersect_STR_df]
	TRF_class_content_df = reduce(merge_dfs, list_merge_TRF_df)
	TRF_class_content_df[['repeat_class', 'repeat_content']] = TRF_class_content_df.apply(calculate_TRF_class_content, axis=1)
	TRF_class_content_df = TRF_class_content_df[columns_out]

	# merge
	sv_repeat_class_df = pd.concat([TEonly_class_content_df, Satelliteonly_class_content_df, TRF_class_content_df], ignore_index=True)
	sv_repeat_class_df['repeat_class'] = sv_repeat_class_df['repeat_class'].astype(int)
	sv_repeat_class_df = sv_repeat_class_df[['sv_id', 'sv_type', 'sv_length', 'repeat_class', 'repeat_content']]

	return sv_repeat_class_df

def sv_class_based_onlyTRF(sv_bed, trf_bed):
	svBed_intersect_TR = (sv_bed.intersect(trf_bed, wao=True).sort()
					   .to_dataframe(names=['chr', 'start', 'end', 'sv_id', 'sv_type', 
							 'TR_chr', 'TR_start', 'TR_end','overlap_length'], 
					  dtype={'chr': str}))
	
	svBed_intersect_TR['sv_length'] = svBed_intersect_TR['end'] - svBed_intersect_TR['start']
	svBed_intersect_TR['TR_length'] = svBed_intersect_TR['TR_end'] - svBed_intersect_TR['TR_start']

	key_list = ['sv_id', 'sv_length']
	svBed_intersect_TR = (svBed_intersect_TR.groupby(key_list)[['overlap_length', 'TR_length']].sum().reset_index())
	svBed_intersect_TR['TR_content'] = svBed_intersect_TR['overlap_length'] / svBed_intersect_TR['sv_length']
	
	columns_out = ['sv_id', 'TR_length', 'TR_content']
	svBed_intersect_TR = svBed_intersect_TR[columns_out]

	return svBed_intersect_TR


def parse_genmap_out(filepath):
	with open(filepath, 'r') as file:
		mappab_dict = {}
		current_chr_name = None
		current_chr_mappab = []

		for line in file:
			line = line.strip()
			if line.startswith('>'):
				if current_chr_name is not None:
					mappab_dict[current_chr_name] = current_chr_mappab
					
				current_chr_name = line[1:]
				current_chr_mappab = []
			else:
				array1 = np.array(line.split(), dtype=np.int64)
				array1[array1>255] = 255
				array1 = np.array(array1, dtype=np.uint8)
				current_chr_mappab = array1

		if current_chr_name is not None:
			mappab_dict[current_chr_name] = current_chr_mappab

	return mappab_dict


def out_mappab(mappab_dict, bed):

	rows_list = []
	for feature in bed:
		chr = str(feature.chrom)
		start = int(feature.start)
		end = int(feature.end)
		id = feature.name
		mappab_median = np.median(mappab_dict[chr][start:end])
		
		row = {'sv_id': id, 'mappab_median': mappab_median}
		rows_list.append(row)

	mappab_df = pd.DataFrame(rows_list, columns=['sv_id', 'mappab_median'])
	
	return mappab_df

def main(args=None):
	parser = argparse.ArgumentParser(description="Extract SV features from different software outputs.")

	parser.add_argument("--ref_sv_vcf", type=str, required=True, help="The VCF file based ref genome", metavar="file")
	parser.add_argument("--alt_sv_bed", type=str, required=True, help="The BED file based alt genome", metavar="file")

	parser.add_argument("--ref_rm", type=str, required=True, help="The ref genome RepeatMasker .out file", metavar='file')
	parser.add_argument("--alt_rm", type=str, required=True, help="The alt genome RepeatMasker .out file", metavar='file')   

	parser.add_argument("--ref_trf", type=str, required=True, help="The ref genome TRF output .gff file", metavar='file')
	parser.add_argument("--alt_trf", type=str, required=True, help="The alt genome TRF output .gff file", metavar='file')
	 
	parser.add_argument("--ref_genmap", type=str, required=True, help="The ref genome genmap output .txt file", metavar='file')
	parser.add_argument("--alt_genmap", type=str, required=True, help="The alt genome genmap output .txt file", metavar='file')

	parser.add_argument("--ref_biser", type=str, required=True, help="The ref genome BISER .out file", metavar='file')
	parser.add_argument("--alt_biser", type=str, required=True, help="The alt genome BISER .out file", metavar='file')

	parser.add_argument("-o", "--out", type=str, default="sv_feature.tsv", help="The output file name, default is: sv_feature.tsv", metavar='file')
	 

	parsed_args = parser.parse_args(args=args)
	ref_vcf_file = parsed_args.ref_sv_vcf
	alt_bed_file = parsed_args.alt_sv_bed

	ref_genmap_file = parsed_args.ref_genmap
	alt_genmap_file = parsed_args.alt_genmap

	ref_RM_file = parsed_args.ref_rm
	alt_RM_file = parsed_args.alt_rm
	 
	ref_biser_file = parsed_args.ref_biser
	alt_biser_file = parsed_args.alt_biser
	 
	ref_trf_file = parsed_args.ref_trf
	alt_trf_file = parsed_args.alt_trf
	out_file = parsed_args.out

	ref_vcf = pysam.VariantFile(ref_vcf_file, "r")
	ref_bed_del = process_vcf_to_del_bed(ref_vcf)
	ref_vcf.close()

	alt_bed_ins = pybedtools.BedTool(alt_bed_file).filter(filter_del).saveas('filtered_alt_bed_file.bed')
	alt_bed_ins = pybedtools.BedTool('filtered_alt_bed_file.bed')
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "1. SVs information has been read.")

	Ref_RM_bed = {}
	Alt_RM_bed = {}

	dtype_spec = {4: str, 5: str, 6: str, 10: str}
	ref_RM = pd.read_csv(ref_RM_file, sep='\s+', header=None, usecols=[4, 5, 6, 10], dtype=dtype_spec)
	alt_RM = pd.read_csv(alt_RM_file, sep='\s+', header=None, usecols=[4, 5, 6, 10], dtype=dtype_spec)

	for repeat_type in ['LINE', 'SINE', 'LTR', 'DNA', 'Satellite']:
		Ref_RM_bed[repeat_type] = process_repeat_masker(ref_RM, repeat_type)
		Alt_RM_bed[repeat_type] = process_repeat_masker(alt_RM, repeat_type)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "2. RepeatMasker out file has been read.")

	ref_SDs_bed = process_biser_out(ref_biser_file, Ref_RM_bed['Satellite'])
	alt_SDs_bed = process_biser_out(alt_biser_file, Alt_RM_bed['Satellite'])
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "3. BISER out file has been read.")

	ref_trf_bed, ref_trf_VNTR_bed, ref_trf_STR_bed = process_TRF_out(ref_trf_file)
	alt_trf_bed, alt_trf_VNTR_bed, alt_trf_STR_bed = process_TRF_out(alt_trf_file)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "4. TRF out file has been read.")

	ref_del_sd_class = sv_class_based_SDs(ref_bed_del, ref_SDs_bed)
	alt_ins_sd_class = sv_class_based_SDs(alt_bed_ins, alt_SDs_bed)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "5. The SDs information has been annotated.")

	ref_del_repeat_class = sv_class_based_Repeats(ref_bed_del, Ref_RM_bed, ref_trf_VNTR_bed, ref_trf_STR_bed)
	alt_ins_repeat_class = sv_class_based_Repeats(alt_bed_ins, Alt_RM_bed, alt_trf_VNTR_bed, alt_trf_STR_bed)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "6. The Repeats information has been annotated.")

	ref_del_TR_class = sv_class_based_onlyTRF(ref_bed_del, ref_trf_bed)
	alt_ins_TR_class = sv_class_based_onlyTRF(alt_bed_ins, alt_trf_bed)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "7. The Tandom repeat information has been annotated.")

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "8. Processing GenMap output file, please be patient...")
	ref_mappab_dict = parse_genmap_out(ref_genmap_file)
	alt_mappab_dict = parse_genmap_out(alt_genmap_file)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "8. GenMap out file has been read.")

	ref_del_mappab  = out_mappab(ref_mappab_dict, ref_bed_del)
	alt_ins_mappab  = out_mappab(alt_mappab_dict, alt_bed_ins)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "9. The Mappability information has been annotated.")

	ref_del_ann_df = pd.merge(pd.merge(pd.merge(ref_del_repeat_class, ref_del_sd_class, on='sv_id', how='inner'),
									ref_del_TR_class, on='sv_id', how='inner'),ref_del_mappab, on='sv_id', how='inner')

	alt_ins_ann_df = pd.merge(pd.merge(pd.merge(alt_ins_repeat_class, alt_ins_sd_class, on='sv_id', how='inner'),
									alt_ins_TR_class, on='sv_id', how='inner'),alt_ins_mappab, on='sv_id', how='inner')

	sv_ann_df = pd.concat([ref_del_ann_df, alt_ins_ann_df], ignore_index=True).round(4)
	sv_ann_df_unique = sv_ann_df.drop_duplicates(subset=['sv_id'], keep='first')
	sv_ann_df_unique.to_csv(out_file, sep='\t', index=False)
	os.remove('filtered_alt_bed_file.bed')

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "10. All annotation tasks have been completed.")
	 

if __name__ == "__main__":
	main()