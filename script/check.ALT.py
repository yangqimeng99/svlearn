import pysam

def check_alt_vcf(alt_sv_vcf, alt_genome, alt_sv_bed):
	vcf = pysam.VariantFile(alt_sv_vcf)
	genome = pysam.FastaFile(alt_genome)

	bed_dict = {}

	with open(alt_sv_bed, 'r') as bed_file:
		for line in bed_file:
			bed_chrom, bed_start, bed_stop, bed_id, bed_svtype = line.strip().split("\t")
			bed_dict[bed_id] =[bed_chrom, int(bed_start), int(bed_stop), bed_svtype]

	for rec in vcf:
		vcf_seq = rec.ref.upper()
		vcf_pos_Seq = genome.fetch(rec.chrom, rec.start, rec.stop).upper()
		bed_pos_Seq = genome.fetch(bed_dict[rec.id][0], bed_dict[rec.id][1], bed_dict[rec.id][2]).upper()

		if vcf_seq != vcf_pos_Seq:
			print("VCF sequence does not match vcf pos genome sequence for SV: " + rec.id)
		elif vcf_seq != bed_pos_Seq:
			print("VCF sequence does not match bed sequence for SV: " + rec.id)
		elif vcf_pos_Seq != bed_pos_Seq:
			print("VCF pos sequence does not match bed sequence for SV: " + rec.id)

	vcf.close()
	genome.close()

alt_sv_vcf = "/svlearn/test_filter/alt_sorted_format_filtered_sv.vcf"
alt_genome = "/github/svlearn/test_filter/alt.fasta"
alt_sv_bed = "/github/svlearn/test_filter/alt_sorted_format_filtered_sv.bed"


check_alt_vcf(alt_sv_vcf, alt_genome, alt_sv_bed)

