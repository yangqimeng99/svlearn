
import os
import re
import sys
import gzip
import pysam
import hashlib
import platform
import argparse
import pandas as pd
from pyfaidx import Fasta
from datetime import datetime

def check_vcf(vcf, ref_genome, output_errors_vcf, output_filtered_vcf):
	sv_ids = set()
	errors = []  # Store tuples of error SVs and records
	filtered_records = []  # Store records without errors SVs
	invalid_base_pattern = re.compile('[^AGCT]')

	for rec in vcf:
		error_messages = []
		svtype = rec.info.get('SVTYPE')
		if isinstance(svtype, tuple):
			svtype = svtype[0]
		# Check if SVTYPE is valid
		if svtype not in ['INS', 'DEL']:
			error_messages.append(f"Invalid_SVTYPE-{svtype}")

		# Check if REF sequence is consistent with the reference genome
		ref_seq = ref_genome.fetch(rec.chrom, rec.pos - 1, rec.stop).upper()
		if rec.ref.upper() != ref_seq:
			error_messages.append("REF_sequence_inconsistent_with_genome_sequence")

		# Check if REF and ALT sequences have different padding bases
		if rec.ref[0].upper() != rec.alts[0][0].upper():
			error_messages.append("Different_padding_base_for_REF_and_ALT_sequence")
		
		# Check if REF and ALT sequences contain invalid bases
		sequences = {'REF': rec.ref.upper(), 'ALT': rec.alts[0].upper()}
		for seq_name, seq in sequences.items():
			if invalid_base_pattern.search(seq):
				error_messages.append(f"Non-ACGT_base_in_{seq_name}_sequence")

		sv_id = rec.id
		if not sv_id or sv_id in sv_ids:
			error_messages.append("SV_ID_either_does_not_exist_or_is_not_unique")
		else:
			sv_ids.add(sv_id)
		
		if error_messages:
			errors.append((rec, error_messages))
		else:
			filtered_records.append(rec)

	new_header = pysam.VariantHeader()
	for rec in vcf.header.records:
		if rec.type == "INFO" and rec["ID"] == "SVTYPE":
			new_header.info.add("SVTYPE", 1, "String", "Type of structural variation")
		else:	
			new_header.add_record(rec)

	for sample in vcf.header.samples:
		new_header.add_sample(sample)

	with pysam.VariantFile(output_filtered_vcf, 'w', header=new_header) as filtered_vcf:
		for rec in filtered_records:
			filtered_vcf.write(rec)
	

	if len(errors) > 0:
		new_header.info.add("ERRORS", 1, "String", "Errors identified in this variant")
		vcf.header.add_line("##INFO=<ID=ERRORS,Number=1,Type=String,Description='Errors identified in this variant'>")
		with pysam.VariantFile(output_errors_vcf, 'w', header=new_header) as errors_vcf:
			for rec, error_msgs in errors:
				rec.info["ERRORS"] = ",".join(error_msgs)
				errors_vcf.write(rec)

	print(f"Finished processing. Errors found in {len(errors)} records.")
	print(f"Output Errors records to: {output_errors_vcf}")
	print(f"Output clean records to: {output_filtered_vcf}")
	print("\n")


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


def generate_contig_name(chr, sv_id):
    # Replace any '-' with '_'
    chr = chr.replace("-", "_")
    sv_id = sv_id.replace("-", "_")
    contig_name = f"{chr}_{sv_id}"

    if len(contig_name) > 40:
        chr_part = chr[:10]
        sv_id_part = sv_id[:10]
        hash_object = hashlib.md5(contig_name.encode())
        short_hash = hash_object.hexdigest()[:8]
        contig_name = f"{chr_part}_{sv_id_part}_{short_hash}"
        # Ensure the contig_name is not longer than 40 characters
        contig_name = contig_name[:40]
    return contig_name

def intervals_overlap(chrom1, start1, end1, chrom2, start2, end2):
	"""Check if two intervals overlap."""
	if chrom1 != chrom2:
		return False
	return start1 <= end2 and start2 <= end1

def process_vcf_overlap(input_vcf):
	"""Process a VCF file to detect overlapping intervals.

	Args:
		input_vcf (str): Path to the input VCF file.

	Returns:
		tuple: Two lists containing non-overlapping and overlapping intervals.
	"""
	# Open the VCF file
	vcf_in = pysam.VariantFile(input_vcf, 'r')
	
	bed_intervals = []          # List to store non-overlapping intervals
	overlapping_intervals = []  # List to store overlapping intervals
	
	for record in vcf_in:
		chrom = record.chrom
		start = record.start
		end = record.stop
	
		current_interval = (chrom, start, end, record)
	
		if not bed_intervals:
			bed_intervals.append(current_interval)
		else:
			last_interval = bed_intervals[-1]
			if intervals_overlap(chrom, start, end, last_interval[0], last_interval[1], last_interval[2]):
				overlapping_intervals.append(current_interval)
			else:
				# If not overlapping, add to bed_intervals
				bed_intervals.append(current_interval)
	
	vcf_in.close()
	return bed_intervals, overlapping_intervals

def extract_alt_genome_from_all_SVs(input_fasta, output_alt_genome, output_alt_bed, bed_intervals, overlapping_intervals, read_length):
	"""Generate an alternative genome by replacing SV regions with alternative sequences.

	Args:
		input_filtered_vcf (str): Path to the input filtered VCF file.
		input_fasta (str): Path to the input reference genome FASTA file.
		output_alt_genome (str): Path to the output alternative genome FASTA file.
		output_alt_bed (str): Path to the output BED file recording SV positions.
	"""
	# Read the reference genome
	ref = Fasta(input_fasta, rebuild=False)
	ref_genome = pysam.FastaFile(input_fasta)
	contigs_names_in_fasta = ref_genome.references
	ref_genome.close()

	# Get the list of contigs with SVs
	contigs_with_sv = set()
	for interval in bed_intervals + overlapping_intervals:
		contigs_with_sv.add(interval[0])

	# Get the list of contigs without SVs
	no_sv_contigs = [seq for seq in contigs_names_in_fasta if seq not in contigs_with_sv]

	out_fa_file = output_alt_genome
	out_bed_file = output_alt_bed

	with open(out_fa_file, "w") as out_fa, open(out_bed_file, "w") as out_bed:
		# Process chromosomes without SVs
		for seq in no_sv_contigs:
			out_fa.write(f">{seq}\n")
			out_fa.write(str(ref[seq]))
			out_fa.write("\n")

		# Process chromosomes with SVs
		for chrom in contigs_names_in_fasta:
			if chrom not in contigs_with_sv:
				continue  # Skip chromosomes without SVs

			# Get non-overlapping SVs for this chromosome
			chrom_bed_intervals = [i for i in bed_intervals if i[0] == chrom]
			if not chrom_bed_intervals:
				# No non-overlapping SVs in this chromosome, write the entire chromosome
				out_fa.write(f">{chrom}\n")
				out_fa.write(str(ref[chrom]))
				out_fa.write("\n")
				continue

			# Sort the SVs by start position
			chrom_bed_intervals.sort(key=lambda x: x[1])

			# Initialize variables
			current_chr = chrom
			seq_len = 0
			end = 0

			# Write chromosome header
			out_fa.write(f">{chrom}\n")

			# Process non-overlapping SVs
			for interval in chrom_bed_intervals:
				chr, start, stop, record = interval
				ref_seq = record.ref
				alt_seq = record.alts[0]
				sv_id = record.id
				svtype = record.info['SVTYPE']

				# Write the sequence before the SV
				sv_left_seq = ref[chr][end:start].seq
				out_fa.write(sv_left_seq)
				seq_len += len(sv_left_seq)
				pseudo_sta = seq_len

				# Write the alternative sequence
				out_fa.write(alt_seq)
				seq_len += len(alt_seq)
				pseudo_end = seq_len

				# Write to alt_bed
				out_bed.write(f"{chr}\t{pseudo_sta}\t{pseudo_end}\t{sv_id}\t{svtype}\n")

				# Update 'end'
				end = start + len(ref_seq)

			# Write the remaining sequence after the last SV
			out_fa.write(ref[chr][end:].seq)
			out_fa.write("\n")

		# Process overlapping SVs
		for interval in overlapping_intervals:
			chr, start, stop, record = interval
			ref_seq = record.ref
			alt_seq = record.alts[0]
			sv_id = record.id
			svtype = record.info['SVTYPE']

			# Extract one read leangth (150 bp) before and after the SV site
			left_pos = max(0, start - read_length)
			right_pos = min(stop + read_length, len(ref[chr]))  # 'stop' is 1-based
			sv_left_seq = ref[chr][left_pos:start].seq
			sv_right_seq = ref[chr][stop:right_pos].seq

			# Create a new contig
			contig_name = generate_contig_name(chr, sv_id)
			out_fa.write(f">{contig_name}\n")

			seq_len = 0

			out_fa.write(sv_left_seq)
			seq_len += len(sv_left_seq)
			pseudo_sta = seq_len

			out_fa.write(alt_seq)
			seq_len += len(alt_seq)
			pseudo_end = seq_len

			out_fa.write(sv_right_seq)
			seq_len += len(sv_right_seq)
			out_fa.write("\n")

			out_bed.write(f"{contig_name}\t{pseudo_sta}\t{pseudo_end}\t{sv_id}\t{svtype}\n")

	ref.close()

def filter_variants(input_vcf, out_clean_vcf, out_close_vcf, slop):
	vcf = pysam.VariantFile(input_vcf)
	variants = []
	for record in vcf:
		if record.info['SVTYPE'] in ['INS', 'DEL']:
			variants.append([record.chrom, record.start, record.stop, record.id])
	
	vcf.close()

	df_variants = pd.DataFrame(variants, columns=['chrom', 'start', 'end', 'id'])

	df_variants.sort_values(by=['chrom', 'start'], inplace=True)
	df_variants['distance_to_next'] = df_variants.groupby('chrom')['start'].shift(-1) - df_variants['end']
	df_variants['distance_to_previous'] = df_variants['start'] - df_variants.groupby('chrom')['end'].shift(1)
	
	filtered_variants = df_variants[((df_variants['distance_to_next'] > slop) | 
									 (df_variants['distance_to_next'].isnull())) & 
									 ((df_variants['distance_to_previous'] > slop) | 
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
	# yes_sv_contigs = [seq for seq in contigs_names_in_fasta if seq in contigs_names_in_vcf]

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
		print(f"    Index file {fai_file} exists or was created successfully.")
	except Exception as e:
		print(f"    Failed to create index for {fasta_file}: {e}")
		sys.exit(1)

def prepare_vcf_header(output_vcf, fasta_file, ref_vcf, all_sv=False):
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
		
		if all_sv:
			vcf_out.write("##INFO=<ID=OVERLAP,Number=1,Type=String,Description='Indicates whether this SV overlaps with other SVs (YES/NO)'>\n")

		vcf_out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def create_alt_vcf(in_ref_vcf, alt_sv_bed, alt_fasta, output_vcf, bed_intervals=None, overlapping_intervals=None):
	vcf_in = pysam.VariantFile(in_ref_vcf, "r")
	
	fasta = pysam.FastaFile(alt_fasta)
	bed_dict = {}

	with open(alt_sv_bed, 'r') as file:
		for line in file:
			parts = line.strip().split()
			chrom, start, end, sv_id, svtype = str(parts[0]), int(parts[1]), int(parts[2]), parts[3], parts[4]
			bed_dict[sv_id] = [chrom, start, end, svtype]

	svtype_mapping = {"INS": "DEL", "DEL": "INS"}

	if overlapping_intervals is not None:
		overlapping_records_set = {item[3].id for item in overlapping_intervals if item[3].id}

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

			if overlapping_intervals is not None:
				overlap = "YES" if record.id in overlapping_records_set else "NO"
				vcf_row = f"{chrom}\t{pos}\t{sv_id}\t{alt_ref_seq}\t{alt_alt_seq}\t{sv_quality}\tPASS\tSVTYPE={svtype};SVLEN={svlen};END={end};OVERLAP={overlap}\n"
			else:
				vcf_row = f"{chrom}\t{pos}\t{sv_id}\t{alt_ref_seq}\t{alt_alt_seq}\t{sv_quality}\tPASS\tSVTYPE={svtype};SVLEN={svlen};END={end}\n"
			
			vcf_out.write(vcf_row)


def get_output_file_path(output_dir_path, file_name):
	return os.path.join(output_dir_path, file_name)

def main(args=None):
	parser = argparse.ArgumentParser(description="from ref genome and vcf create alt genome and vcf.")

	parser.add_argument("--ref_fasta", type=str, required=True, help="The Fasta file of ref genome", metavar="file")
	parser.add_argument("--ref_sv_vcf", type=str, required=True, help="The SV VCF file based ref genome", metavar="file")
	parser.add_argument('--no-filter-overlaps', action='store_true', help='If specified, use all SVs to construct the ALT genome; otherwise, filter out overlapping SVs to improve genotyping performance.')
	parser.add_argument("--min_distance", type=str, default="300", help="Filter out pairs of SVs that are less than this distance apart, default: 300", metavar="int")
	parser.add_argument("--read_length", type=str, default="150", help="Short reads leagth, default: 150", metavar="int")
	parser.add_argument("-o", "--out", type=str, default="prepareAlt_output", help="The output dir name of prepareAlt", metavar="dir")


	parsed_args = parser.parse_args(args=args)
	input_fasta = parsed_args.ref_fasta
	input_vcf = parsed_args.ref_sv_vcf
	output_dir = parsed_args.out
	slop = max(int(parsed_args.min_distance), 0)
	read_length = int(parsed_args.read_length)

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "python version:", platform.python_version())
	print(current_time, "short reads length:", read_length)

	output_dir_path = os.path.abspath(output_dir)
	os.makedirs(output_dir_path, exist_ok=True)

	vcf = pysam.VariantFile(input_vcf)
	fasta = pysam.FastaFile(input_fasta)

	output_errors_vcf = get_output_file_path(output_dir_path, "ref_errors.vcf")
	output_clean_vcf = get_output_file_path(output_dir_path, "ref_clean.vcf")
	check_vcf(vcf, fasta, output_errors_vcf, output_clean_vcf)
	vcf.close()
	fasta.close()

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "00. All records meet the conditions, the program continues to run.")

	output_sorted_vcf_path = get_output_file_path(output_dir_path, "ref_sorted.vcf")
	sort_vcf(output_clean_vcf, output_sorted_vcf_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "01. VCF sorted successfully!")
	print("    The ref sorted vcf is:", output_sorted_vcf_path)

	output_format_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format.vcf")
	process_vcf(output_sorted_vcf_path, output_format_vcf_path)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "02. VCF formatted successfully!")
	print("    The ref sorted format vcf is:", output_format_vcf_path)

	if parsed_args.no_filter_overlaps:
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "03. use all SVs to construct the ALT genome")
		output_filtered_vcf_path = output_format_vcf_path
		output_alt_genome_path = get_output_file_path(output_dir_path, "alt.fasta")   # Output path for the alternative genome FASTA file
		output_alt_bed_path = get_output_file_path(output_dir_path, "alt.bed")     # Output path for the BED file
		
		# Process the VCF file to get non-overlapping and overlapping intervals
		bed_intervals, overlapping_intervals = process_vcf_overlap(output_format_vcf_path)
		extract_alt_genome_from_all_SVs(input_fasta, output_alt_genome_path, output_alt_bed_path, bed_intervals, overlapping_intervals, read_length)
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "04. ALT genome created successfully!")
		print("    The ALT genome Fasta is:", output_alt_genome_path)
		print("    The ALT genome SVs BED is:", output_alt_bed_path)

		out_alt_vcf_file_path = get_output_file_path(output_dir_path, "alt_sorted_format.vcf")
		create_fasta_index(output_alt_genome_path)
		prepare_vcf_header(out_alt_vcf_file_path, output_alt_genome_path, output_filtered_vcf_path, all_sv=True)
		create_alt_vcf(output_filtered_vcf_path, output_alt_bed_path, output_alt_genome_path, out_alt_vcf_file_path, bed_intervals, overlapping_intervals)

	else:
		output_filtered_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format_filtered_sv.vcf")
		output_close_vcf_path = get_output_file_path(output_dir_path, "ref_sorted_format_close_sv.vcf")
		filter_variants(output_format_vcf_path, output_filtered_vcf_path, output_close_vcf_path, slop)
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
