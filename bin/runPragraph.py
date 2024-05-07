# %%
import re
import os
import sys
import glob
import math
import pysam
import shutil
import argparse
import pybedtools
import subprocess
import numpy as np
import polars as pl
from pyfaidx import Fasta
from functools import partial
from datetime import datetime

# %%
def get_output_file_path(output_dir_path, file_name):
    return os.path.join(output_dir_path, file_name)

# %%
def main(args=None):
    parser = argparse.ArgumentParser(description="Run paragraph in ref and alt bam.")
    parser.add_argument("--ref_fasta", type=str, required=True, help="The ref genome Fasta file", metavar="file")
    parser.add_argument("--alt_fasta", type=str, required=True, help="The alt genome Fasta file", metavar="file")
    parser.add_argument("--ref_sv_vcf", type=str, required=True, help="The VCF file based ref genome", metavar="file")
    parser.add_argument("--alt_sv_vcf", type=str, required=True, help="The BED file based alt genome", metavar="file")
    parser.add_argument("--ref_bam", type=str, required=True, help="The ref genome bam file", metavar="file")
    parser.add_argument("--alt_bam", type=str, required=True, help="The alt genome bam file", metavar="file")
    parser.add_argument("-t", "--threads", type=str, default="1", help="threads", metavar="int")
    parser.add_argument("-o", "--out", type=str, default="Para_feature_out", help="The output dir name, Default: Para_feature_out", metavar='file')

    parsed_args = parser.parse_args(args=args)
    ref_vcf_file = parsed_args.ref_sv_vcf
    alt_vcf_file = parsed_args.alt_sv_vcf
    ref_bam_file = parsed_args.ref_bam
    alt_bam_file = parsed_args.alt_bam
    ref_fasta_file = parsed_args.ref_fasta
    alt_fasta_file = parsed_args.alt_fasta
    threads = parsed_args.threads
    output_dir = parsed_args.out

    output_dir_path = os.path.abspath(output_dir)
    os.makedirs(output_dir_path, exist_ok=True)
    ref_idxdepth_json = get_output_file_path(output_dir_path, "ref.idxdepth.json")
    alt_idxdepth_json = get_output_file_path(output_dir_path, "alt.idxdepth.json")

    # Run idxdepth
    script_dir = os.path.dirname(os.path.realpath(__file__))
    idxdepth_path = os.path.join(script_dir, 'paragraph-v2.4a/bin/idxdepth')

    idxdepth_ref_command = [
        idxdepth_path,
        '-b', ref_bam_file,
        '-r', ref_fasta_file,
        '-o', ref_idxdepth_json,
        '--threads', threads,
        '--log-file', f"{ref_idxdepth_json}.idxdepth.log"]

    try:
        subprocess.run(idxdepth_ref_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running idxdepth in ref bam: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred while running idxdepth in ref bam: {e}")
        sys.exit(1)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "01. Run idxdepth in ref bam done.")

    idxdepth_alt_command = [
        idxdepth_path,
        '-b', alt_bam_file,
        '-r', alt_fasta_file,
        '-o', alt_idxdepth_json,
        '--threads', threads,
        '--log-file', f"{alt_idxdepth_json}.log"]

    try:
        subprocess.run(idxdepth_alt_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running idxdepth in alt bam: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred while running idxdepth in alt bam: {e}")
        sys.exit(1)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "02. Run idxdepth in alt bam done.")

    # make manifest file
    ref_idxdepth_df = pl.read_json(ref_idxdepth_json)
    alt_idxdepth_df = pl.read_json(alt_idxdepth_json)
    ref_depth = float(ref_idxdepth_df["autosome"][0]["depth"])
    alt_depth = float(alt_idxdepth_df["autosome"][0]["depth"])
    ref_read_len = int(ref_idxdepth_df["read_length"][0])
    alt_read_len = int(alt_idxdepth_df["read_length"][0])
    ref_id = "GT_paragraph_ref"
    alt_id = "GT_paragraph_alt"
    ref_manifest = get_output_file_path(output_dir_path, "ref.manifest")
    alt_manifest = get_output_file_path(output_dir_path, "alt.manifest")

    with open(ref_manifest, "w") as ref_manifest_out, open(alt_manifest, "w") as alt_manifest_out:
        ref_manifest_out.write(f"id\tpath\tdepth\tread length\n")
        ref_manifest_out.write(f"{ref_id}\t{ref_bam_file}\t{ref_depth}\t{ref_read_len}\n")
        alt_manifest_out.write(f"id\tpath\tdepth\tread length\n")
        alt_manifest_out.write(f"{alt_id}\t{alt_bam_file}\t{alt_depth}\t{alt_read_len}\n")

    # Run paragraph
    paragraph_path = os.path.join(script_dir, 'paragraph-v2.4a/bin/multigrmpy.py')
    ref_max_depth = int(20 * ref_depth / 1)
    alt_max_depth = int(20 * alt_depth / 1)
    ref_paragraph_out = get_output_file_path(output_dir_path, "ref_paragraph_out")
    alt_paragraph_out = get_output_file_path(output_dir_path, "alt_paragraph_out")
    ref_paragraph_out_tmp = get_output_file_path(ref_paragraph_out, "tmp")
    alt_paragraph_out_tmp = get_output_file_path(alt_paragraph_out, "tmp")

    ref_paragraph_command = [
        sys.executable, paragraph_path,
        '-i', ref_vcf_file,
        '-M', str(ref_max_depth),
        '-m', ref_manifest,
        '-r', ref_fasta_file,
        '-t', threads,
        '-o', ref_paragraph_out,
        '--scratch-dir', ref_paragraph_out_tmp
    ]

    try:
        subprocess.run(ref_paragraph_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running paragraph in ref bam: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred while running paragraph in ref bam: {e}")
        sys.exit(1)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "03. Run paragraph in ref bam done.")


    alt_paragraph_command = [
        sys.executable, paragraph_path,
        '-i', alt_vcf_file,
        '-M', str(alt_max_depth),
        '-m', alt_manifest,
        '-r', alt_fasta_file,
        '-t', threads,
        '-o', alt_paragraph_out,
        '--scratch-dir', alt_paragraph_out_tmp
    ]

    try:
        subprocess.run(alt_paragraph_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running paragraph in alt bam: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred while running paragraph in alt bam: {e}")
        sys.exit(1)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "04. Run paragraph in alt bam done.")

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

    ref_paragraph_vcf_file = get_output_file_path(ref_paragraph_out, "genotypes.vcf.gz")
    alt_paragraph_vcf_file = get_output_file_path(alt_paragraph_out, "genotypes.vcf.gz")
    
    # get ref info
    ref_paragraph_gt_info = []
    ref_paragraph_vcf = pysam.VariantFile(ref_paragraph_vcf_file)
    for record in ref_paragraph_vcf:
        sv_id = record.id
        ref_para_geno_info = record.samples['GT_paragraph_ref']
        ref_sv_gt_tup = ref_para_geno_info['GT']
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

    # get alt info
    alt_paragraph_gt_info = []
    alt_paragraph_vcf = pysam.VariantFile(alt_paragraph_vcf_file)
    for record in alt_paragraph_vcf:
        sv_id = record.id
        alt_para_geno_info = record.samples['GT_paragraph_alt']
        alt_sv_gt_tup = alt_para_geno_info['GT']
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
    alt_paragraph_vcf.close()
    para_info_df = ref_info_df.join(alt_info_df, on="sv_id", how="left")

    para_info_df = (para_info_df
				.with_columns((pl.col("ref_DP")-pl.col("alt_DP")).alias("DP"))
				.with_columns((pl.col("ref_PL")-pl.col("alt_PL")).alias("PL")))

    para_info_out = get_output_file_path(output_dir_path, "para_feature.tsv")
    para_info_df.select(['sv_id', 'ref_GT', 'alt_GT', 'ref_FT', 'alt_FT', 'DP', 'PL']).write_csv(para_info_out, separator="\t")
    shutil.rmtree(ref_paragraph_out_tmp)
    shutil.rmtree(alt_paragraph_out_tmp)
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "05. The tmp directory generated by running paragraph has been deleted.")
    print("All done.")

if __name__ == "__main__":
    main()
