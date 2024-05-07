
import os
import sys
import gzip
import pysam
import platform
import argparse
import pandas as pd
from pyfaidx import Fasta
from datetime import datetime

def check_vcf(vcf, ref_genome):
	sv_ids = set()
	for rec in vcf:
		svtype = rec.info.get('SVTYPE')
		if svtype not in ['INS', 'DEL']:
			raise ValueError(f"Error: At position {rec.chrom}:{rec.pos}, SV type is not INS or DEL.")

		ref_seq = ref_genome.fetch(rec.chrom, rec.pos - 1, rec.stop).upper()
		if rec.ref.upper() != ref_seq:
			raise ValueError(f"Error: At position {rec.chrom}:{rec.pos}, REF sequence is incomplete or inconsistent with the genome sequence.")
		elif rec.ref[0].upper() != rec.alts[0][0].upper():
			raise ValueError(f"Error: At position {rec.chrom}:{rec.pos}, Different padding base for REF and ALT sequence.")

		sv_id = rec.id
		if not sv_id or sv_id in sv_ids:
			raise ValueError(f"Error: At position {rec.chrom}:{rec.pos}, SV ID either does not exist or is not unique.")
		sv_ids.add(sv_id)

def sort_vcf(input_vcf, output_sorted_vcf):
    """
    Sort a VCF file using pysam.
    """
    with pysam.VariantFile(input_vcf) as vcf_in:
        variants = [variant for variant in vcf_in]
        variants.sort(key=lambda x: (x.chrom, x.pos))

        with pysam.VariantFile(output_sorted_vcf, 'w', header=vcf_in.header) as vcf_out:
            for variant in variants:
                vcf_out.write(variant)

def process_vcf(input_vcf_path, output_vcf_path):
	if input_vcf_path.endswith('.gz'):
		opener = gzip.open
	else:
		opener = open

	with opener(input_vcf_path, 'r') as input_vcf, open(output_vcf_path, 'w') as output_vcf:
		for line in input_vcf:
			if line.startswith('##FORMAT'):
				continue
			
			if line.startswith('##'):
				output_vcf.write(line)
				continue

			output_vcf.write('\t'.join(line.strip().split('\t', 8)[:8]) + '\n')

def filter_variants(input_vcf, out_clean_vcf, out_close_vcf, slop=75):
	slop = int(slop)
	vcf = pysam.VariantFile(input_vcf)
	variants = []
	for record in vcf:
		if record.info['SVTYPE'] in ['INS', 'DEL']:
			variants.append([record.chrom, record.start, record.stop, record.id])
	
	vcf.close()

	df_variants = pd.DataFrame(variants, columns=['chrom', 'start', 'end', 'id'])
	df_variants['start_slop'] = df_variants['start'] - slop
	df_variants['end_slop'] = df_variants['end'] + slop
	
	df_variants.sort_values(by=['chrom', 'start_slop'], inplace=True)
	df_variants['distance_to_next'] = df_variants.groupby('chrom')['start_slop'].shift(-1) - df_variants['end_slop']
	df_variants['distance_to_previous'] = df_variants['start_slop'] - df_variants.groupby('chrom')['end_slop'].shift(1)
	
	filtered_variants = df_variants[((df_variants['distance_to_next'] > 150) | 
									 (df_variants['distance_to_next'].isnull())) & 
									 ((df_variants['distance_to_previous'] > 150) | 
									  (df_variants['distance_to_previous'].isnull()))]
	filtered_variants_id_list = filtered_variants['id'].tolist()

	with pysam.VariantFile(input_vcf) as vcf_in:
		with pysam.VariantFile(out_clean_vcf, 'w', header=vcf_in.header) as vcf_clean_out:
			with pysam.VariantFile(out_close_vcf, 'w', header=vcf_in.header) as vcf_close_out:
				for record in vcf_in:
					if record.id in filtered_variants_id_list:
						vcf_clean_out.write(record)
					else:
						vcf_close_out.write(record)

	vcf_in.close() 

def extract_alt_genome(input_filtered_vcf, input_fasta, output_alt_genome, output_alt_bed):

	ref_genome = pysam.FastaFile(input_fasta)
	contigs_names_in_fasta = ref_genome.references
	ref_genome.close()

	ref_vcf = pysam.VariantFile(input_filtered_vcf)
	contigs_names_in_vcf = []
	for rec in ref_vcf:
		contigs_names_in_vcf.append(rec.chrom)

	ref_vcf.close()

	no_sv_contigs = [seq for seq in contigs_names_in_fasta if seq not in contigs_names_in_vcf]
	yes_sv_contigs = [seq for seq in contigs_names_in_fasta if seq in contigs_names_in_vcf]

	current_chr = None
	seq_len = 0
	end = 0

	out_fa_file = output_alt_genome
	out_bed_file = output_alt_bed

	ref = Fasta(input_fasta, rebuild=False)
	first_record_chr = contigs_names_in_vcf[0]
	clean_sv_vcf = pysam.VariantFile(input_filtered_vcf)

	with open(out_fa_file, "w") as out_fa, open(out_bed_file, "w") as out_bed:
		out_fa.write(f">{first_record_chr}\n")
		for record in clean_sv_vcf:
			chr = record.chrom
			start = record.start
			ref_seq = record.ref
			alt_seq = record.alts[0]
			sv_id = record.id
			svtype = record.info['SVTYPE']
			if current_chr == None:
				current_chr = chr
			elif chr != current_chr:
				out_fa.write(ref[current_chr][end:].seq)
				out_fa.write(f"\n>{chr}\n")
				current_chr = chr
				seq_len = 0
				end = 0

			sv_left_seq = ref[chr][end:start].seq
			out_fa.write(sv_left_seq)
			seq_len += len(sv_left_seq)
			pseudo_sta = seq_len

			out_fa.write(alt_seq)
			seq_len += len(alt_seq)
			pseudo_end = seq_len
			out_bed.write(f"{chr}\t{pseudo_sta}\t{pseudo_end}\t{sv_id}\t{svtype}\n")
			end = start + len(ref_seq)

		out_fa.write(ref[current_chr][end:].seq)
		if len(no_sv_contigs) != 0:
			for seq in no_sv_contigs:
				out_fa.write(f"\n>{seq}\n")
				out_fa.write(str(ref[seq]))
		
		out_fa.write("\n")

	ref.close()
	clean_sv_vcf.close()


def create_fasta_index(fasta_file):
	"""Check if FASTA index exists, create one if it doesn't."""
	fai_file = fasta_file + ".fai"
	try:
		Fasta(fasta_file)
		print(f"05. Index file {fai_file} exists or was created successfully.")
	except Exception as e:
		print(f"Failed to create index for {fasta_file}: {e}")
		sys.exit(1)

def prepare_vcf_header(output_vcf, fasta_file, ref_vcf):
	"""Prepare the header of the output VCF file based on the given FASTA and reference VCF."""
	with open(output_vcf, 'w') as vcf_out:
		current_time = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
		vcf_out.write("##fileformat=VCFv4.2\n")
		vcf_out.write(f'##fileDate="{current_time}"\n')
		
		fasta = Fasta(fasta_file)
		for seq in fasta:
			vcf_out.write(f"##contig=<ID={seq.name},length={len(seq)}>\n")

		with open(ref_vcf) as ref_vcf_in:
			for line in ref_vcf_in:
				if line.startswith("##ALT=") or line.startswith("##FILTER=") or line.startswith("##INFO="):
					vcf_out.write(line)

		vcf_out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def create_alt_vcf(in_ref_vcf, alt_sv_bed, alt_fasta, output_vcf):
	vcf_in = pysam.VariantFile(in_ref_vcf, "r")
	
	fasta = pysam.FastaFile(alt_fasta)
	bed_dict = {}

	with open(alt_sv_bed, 'r') as file:
		for line in file:
			parts = line.strip().split()
			chrom, start, end, sv_id, svtype = str(parts[0]), int(parts[1]), int(parts[2]), parts[3], parts[4]
			bed_dict[sv_id] = [chrom, start, end, svtype]

	svtype_mapping = {"INS": "DEL", "DEL": "INS"}

	with open(output_vcf, 'a') as vcf_out:
		for record in vcf_in:
			sv_id = record.id
			chrom, start, end, svtype = bed_dict[sv_id]
			if record.qual is None:
				sv_quality = 60
			else:
				sv_quality = int(record.qual)
			alt_ref_seq = record.alts[0]
			alt_alt_seq = record.ref
			
			chrom, start, end, svtype = bed_dict[sv_id]
			pos = start + 1
			svtype = svtype_mapping[svtype]

			svlen = len(alt_alt_seq) if svtype == "INS" else len(alt_ref_seq)
			
			vcf_row = f"{chrom}\t{pos}\t{sv_id}\t{alt_ref_seq}\t{alt_alt_seq}\t{sv_quality}\tPASS\tSVTYPE={svtype};SVLEN={svlen};END={end}\n"
			vcf_out.write(vcf_row)


def get_output_file_path(output_dir_path, file_name):
	return os.path.join(output_dir_path, file_name)

def main(args=None):
	parser = argparse.ArgumentParser(description="from ref genome and vcf create alt genome and vcf.")

	parser.add_argument("--ref_fasta", type=str, required=True, help="The Fasta file of ref genome", metavar="file")
	parser.add_argument("--ref_sv_vcf", type=str, required=True, help="The SV VCF file based ref genome", metavar="file")
	parser.add_argument("-o", "--out", type=str, default="prepareAlt_output", help="The output dir name of prepareAlt", metavar="dir")

	parsed_args = parser.parse_args(args=args)
	input_fasta = parsed_args.ref_fasta
	input_vcf = parsed_args.ref_sv_vcf
	output_dir = parsed_args.out

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "python version:", platform.python_version())

	vcf = pysam.VariantFile(input_vcf)
	fasta = pysam.FastaFile(input_fasta)
	try:
		check_vcf(vcf, fasta)
	except ValueError as e:
		print("VCF format wrong:", e)
		exit(1)

	vcf.close()
	fasta.close()

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "00. All records meet the conditions, the program continues to run.")

	output_dir_path = os.path.abspath(output_dir)
	os.makedirs(output_dir_path, exist_ok=True)

	output_sorted_vcf_path = get_output_file_path(output_dir_path, "ref_sorted.vcf")
	sort_vcf(input_vcf, output_sorted_vcf_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "01. VCF sorted successfully!")
	print("    The ref sorted vcf is:", output_sorted_vcf_path)

	output_format_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format.vcf")
	process_vcf(output_sorted_vcf_path, output_format_vcf_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "02. VCF formatted successfully!")
	print("    The ref sorted format vcf is:", output_format_vcf_path)

	output_filtered_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format_filtered_sv.vcf")
	output_close_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format_close_sv.vcf")
	filter_variants(output_format_vcf_path, output_filtered_vcf_path, output_close_vcf_path, slop=75)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "03. VCF SVs filtered successfully!")
	print(f"    Filtered variants are written to {output_filtered_vcf_path}")
	print(f"    Close variants are written to {output_close_vcf_path}")

	output_alt_genome_path = get_output_file_path(output_dir_path, "alt.fasta")
	output_alt_bed_path = get_output_file_path(output_dir_path, "alt_sorted_format_filtered_sv.bed")
	extract_alt_genome(output_filtered_vcf_path, input_fasta, output_alt_genome_path, output_alt_bed_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "04. ALT genome created successfully!")
	print("    The ALT genome Fasta is:", output_alt_genome_path)
	print("    The ALT genome SVs BED is:", output_alt_bed_path)

	out_alt_vcf_file_path = get_output_file_path(output_dir_path, "alt_sorted_format_filtered_sv.vcf")
	create_fasta_index(output_alt_genome_path)
	prepare_vcf_header(out_alt_vcf_file_path, output_alt_genome_path, output_filtered_vcf_path)
	create_alt_vcf(output_filtered_vcf_path, output_alt_bed_path, output_alt_genome_path, out_alt_vcf_file_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "05. ALT VCF created successfully!")
	print("    The ALT genome SVs VCF is:", out_alt_vcf_file_path)
	print("ALL DONE!")
	
if __name__ == "__main__":
	main()
