import pysam
import argparse
import pandas as pd

def collect_reads_info(bam, bed_file, bam_allele, reads_dict):
    with open(bed_file, "r") as bed:
        for line in bed:
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            sv_id = line[3]
            sv_type = line[4]
            for read in bam.fetch(chr, start, end):
                if read.is_secondary or read.is_unmapped or read.is_duplicate:
                    continue

                read_name = read.query_name
                readsuffix = '_1' if read.is_read1 else '_2'
                read_id = f'{read_name}{readsuffix}'
                mapping_quality = read.mapping_quality
                start_position = read.reference_start
                end_position = read.reference_end

                if read_id in reads_dict:
                    if mapping_quality > reads_dict[read_id]['mapping_quality']:
                        reads_dict[read_id] = {
                            'start_position': start_position,
                            'end_position': end_position,
                            'mapping_quality': mapping_quality,
                            'bam_allele': bam_allele,
                            'sv_id': sv_id,
                            'sv_type': sv_type
                        }
                else:
                    reads_dict[read_id] = {
                        'start_position': start_position,
                        'end_position': end_position,
                        'mapping_quality': mapping_quality,
                        'bam_allele': bam_allele,
                        'sv_id': sv_id,
                        'sv_type': sv_type
                    }
    return reads_dict

def main(args=None):
    parser = argparse.ArgumentParser(description="collect sv site reads mapping info from ref and alt bam.")
    parser.add_argument("--ref_bam", type=str, required=True, help="The ref genome bam file.")
    parser.add_argument("--alt_bam", type=str, required=True, help="The alt genome bam file.")
    parser.add_argument("--ref_sv_bed", type=str, required=True, help="The ref sv bed file.")
    parser.add_argument("--alt_sv_bed", type=str, required=True, help="The alt sv bed file.")
    parser.add_argument("--prefix", type=str, default="out", help="The output file prefix.")
    
    parsed_args = parser.parse_args(args=args)
    ref_bam_file = parsed_args.ref_bam
    alt_bam_file = parsed_args.alt_bam
    ref_sv_bed_file = parsed_args.ref_sv_bed
    alt_sv_bed_file = parsed_args.alt_sv_bed
    prefix = parsed_args.prefix

    summary_file = f'{prefix}_summary.txt'

    reads_dict = {}
    ref_bam = pysam.AlignmentFile(ref_bam_file, "rb")
    reads_dict = collect_reads_info(ref_bam, ref_sv_bed_file, 0, reads_dict)
    ref_bam.close()
    ref_dict = dict(reads_dict)

    alt_bam = pysam.AlignmentFile(alt_bam_file, "rb")
    reads_dict = collect_reads_info(alt_bam, alt_sv_bed_file, 1, reads_dict)
    alt_bam.close()
    alt_dict = dict(reads_dict)

    df_ref = pd.DataFrame.from_dict(ref_dict, orient='index')
    df_ref = df_ref.reset_index()
    df_ref = df_ref.rename(columns={'index': 'read_id'})

    print(f"the reads number of sv site in ref bam: {len(df_ref)}")
    ref_output_file = f'{prefix}_ref_sv_reads_info.tsv'
    df_ref.to_csv(ref_output_file, sep='\t', index=False)

    df_alt = pd.DataFrame.from_dict(alt_dict, orient='index')
    df_alt = df_alt.reset_index()
    df_alt = df_alt.rename(columns={'index': 'read_id'})

    print(f"the reads number of sv site in ref and alt bam: {len(df_alt)}")
    print(f"the add mapped reads number of sv site: {len(df_alt) - len(df_ref)}")
    alt_output_file = f'{prefix}_ref_alt_sv_reads_info.tsv'
    df_alt.to_csv(alt_output_file, sep='\t', index=False)

    df_add = df_alt["bam_allele"] == 1
    add_best_mapped_reads_out = f'{prefix}_add_best_mapped_sv_reads_info.tsv'
    df_add_best_mapped = df_alt[df_add]
    df_add_best_mapped.to_csv(add_best_mapped_reads_out, sep='\t', index=False)
    print(f"the add best mapped reads number of sv site: {len(df_add_best_mapped)}")

    
    with open(summary_file, "w") as f:
        f.write(f"the reads number of sv site in ref bam:\t{len(df_ref)}\n")
        f.write(f"the reads number of sv site in ref and alt bam:\t{len(df_alt)}\n")
        f.write(f"the add mapped reads number of sv site:\t{len(df_alt) - len(df_ref)}\n")
        f.write(f"the add best mapped reads number of sv site:\t{len(df_add_best_mapped)}\n")

if __name__ == "__main__":
    main()


