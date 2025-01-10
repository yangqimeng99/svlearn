import os
import sys
import gzip
import argparse
import numpy as np
import polars as pl
from functools import reduce
from datetime import datetime

from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from joblib import load

def pre_data(pl_df, change_dtype_dict, check_null, check_nan):

    pl_df = pl_df.with_columns(
        pl.when(pl_df['sv_type'] == "INS").then(0)
        .when(pl_df['sv_type'] == "DEL").then(1)
        .otherwise(pl.lit(np.nan)).alias("sv_type")
        )
    
    # condition1 = reduce(lambda a, b: a | b, [pl_df[col] == "-" for col in check_null])
    condition1 = reduce(lambda a, b: a | b, [pl_df[col].is_null() for col in check_null])
    pl_df_filted = pl_df.filter(~condition1)

    condition2 = reduce(lambda a, b: a | b, [pl_df_filted[col].is_nan() for col in check_nan])
    pl_df_filted = pl_df_filted.filter(~condition2).cast(change_dtype_dict)

    return pl_df_filted, len(pl_df), len(pl_df_filted)

def auto_open(filename, mode='rt'):
    if filename.endswith('.vcf'):
        return open(filename, mode)
    elif filename.endswith('.vcf.gz'):
        if "w" in mode:
            sys.exit('Error: The output File extension must be .vcf')
        else:
            return gzip.open(filename, mode)
    else:
         sys.exit('Error: The input File extension must be vcf or vcf.gz, the output File extension must be .vcf')

def main(args=None):
    parser = argparse.ArgumentParser(description="SV genotyping was performed using a machine learning model.")
    parser.add_argument("-v", "--ref_sv_vcf", type=str, required=True, help="The VCF file based ref genome", metavar="file")
    parser.add_argument("-m", "--model", type=str, required=True, help="The machine learning model file", metavar="file")
    parser.add_argument("-s", "--sv_feature", type=str, required=True, help="SV feature from ref and alt genome, the output of SVfeature", metavar="file")
    parser.add_argument("-a", "--align_feature", type=str, required=True, help="BreakPoint and ReadDepth feature from 2bam, the output of alignFeature", metavar="file")
    parser.add_argument("-p", "--paragraph_feature", type=str, default=None, help="Paragraph feature from 2bam, the output of runParagraph", metavar="file")
    parser.add_argument("-n", "--name", type=str, required=True, help="The sample name in out vcf file", metavar="str")
    parser.add_argument("-o", "--out", type=str, default="svlearn_genotype.vcf", help="The out vcf file, Default: svlearn_genotype.vcf", metavar="file")
    
    parsed_args = parser.parse_args(args=args)
    ref_sv_vcf_file = parsed_args.ref_sv_vcf
    model_file = parsed_args.model
    sv_feature_file = parsed_args.sv_feature
    align_feature_file = parsed_args.align_feature
    paragraph_feature_file = parsed_args.paragraph_feature
    out_name = parsed_args.name
    out_file = parsed_args.out

    sv_feature = pl.read_csv(sv_feature_file, separator='\t')
    align_feature = pl.read_csv(align_feature_file, separator='\t', null_values=["-"])
    paragraph_feature_file = None if paragraph_feature_file == "None" else paragraph_feature_file
    if paragraph_feature_file is not None:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(current_time, "paragraph feature input")
        feature_dtype_change_dict = {
            "sv_type":pl.Int32,
            "BP_lp_homref":pl.Float64, 
            "BP_lp_het":pl.Float64, 
            "BP_lp_homalt":pl.Float64,
            "ref_GT":pl.Int32,
            "alt_GT":pl.Int32,
            "PL":pl.Int32,
            "ref_FT":pl.Int32,
            "alt_FT":pl.Int32
            }
        check_null = ["BP_lp_homref", "BP_lp_het", "BP_lp_homalt"]
        check_nan = ["sv_type", "ref_GT", "alt_GT", "PL"]
        categorical_features = ['sv_type', 'repeat_class', 'SDs_class', 'ref_GT', 'alt_GT', 'ref_FT', 'alt_FT']
        categories_per_feature = {
			'sv_type': [0, 1],
			'repeat_class': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
			'SDs_class': [0, 1],
			'ref_GT': [0, 1, 2],
			'alt_GT': [0, 1, 2],
			'ref_FT': [0, 1, 2, 3, 4, 5, 6],
			'alt_FT': [0, 1, 2, 3, 4, 5, 6]
		}
        paragraph_feature = pl.read_csv(paragraph_feature_file, separator='\t')
        feature_df = sv_feature.join(align_feature, on="sv_id", how="inner").join(paragraph_feature, on="sv_id", how="inner")
    else:
        feature_dtype_change_dict = {
            "sv_type":pl.Int32,
            "BP_lp_homref":pl.Float64, 
            "BP_lp_het":pl.Float64, 
            "BP_lp_homalt":pl.Float64
            }
        check_null = ["BP_lp_homref", "BP_lp_het", "BP_lp_homalt"]
        check_nan = ["sv_type"]
        categorical_features = ['sv_type', 'repeat_class', 'SDs_class']
        categories_per_feature = {
			'sv_type': [0, 1],
			'repeat_class': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
			'SDs_class': [0, 1]
		}
        feature_df = sv_feature.join(align_feature, on="sv_id", how="inner")

    filted_id_df, input_num, filtered_num, = pre_data(feature_df, feature_dtype_change_dict, check_null, check_nan)
    filted_df = filted_id_df.drop("sv_id")
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "The number of SV entries entered is:", input_num)
    print(current_time, "The number of SV entries after filtering is:", filtered_num)

    one_hot_encoders = {
        feature: OneHotEncoder(categories=[categories_per_feature[feature]], handle_unknown='ignore')
        for feature in categorical_features
    }
    transformers = [(feature, one_hot_encoders[feature], [feature]) for feature in categorical_features]
    column_transformer = ColumnTransformer(transformers=transformers, remainder='passthrough')

    model = load(model_file)
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "The model has been loaded successfully")
    
    filted_df_pd = filted_df.to_pandas()
    geno_x_encoded = column_transformer.fit_transform(filted_df_pd)
    pred_y_val = model.predict(geno_x_encoded)

    filted_id_df = (filted_id_df
                    .with_columns(pl.Series("GT_pred", pred_y_val))
                    .select(["sv_id", "GT_pred"]))

    out_df = (feature_df
                    .join(filted_id_df, on='sv_id', how='left', validate='1:1')
                    .select(["sv_id", "GT_pred"])
                    .with_columns(
                        pl.when(pl.col('GT_pred') == 0).then(pl.lit("0/0"))
                        .when(pl.col('GT_pred') == 1).then(pl.lit("0/1"))
                        .when(pl.col('GT_pred') == 2).then(pl.lit("1/1"))
                        .otherwise(pl.lit("./.")).alias("GT_pred")))


    out_dict = out_df.to_dict(as_series=False)
    out_dict = {out_dict["sv_id"][i]: out_dict["GT_pred"][i] for i in range(len(out_df))}
    vcf_out_file = out_file

    with auto_open(vcf_out_file, "w") as vcf_out:
        with auto_open(ref_sv_vcf_file, "rt") as ref_vcf:
            vcf_out.write("##fileformat=VCFv4.2\n")
            for line in ref_vcf:
                if line.startswith("##ALT=") or line.startswith("##FILTER=") or line.startswith("##INFO=") or line.startswith("##contig="):
                    vcf_out.write(line)

            vcf_out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf_out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{out_name}\n")
                  
        with auto_open(ref_sv_vcf_file, "rt") as ref_vcf:
            for line in ref_vcf:
                if line[0] != "#":
                    items = line.strip().split()
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = items[:9]
                    gt_pred = out_dict.get(ID, "./.")
                    outitems = [CHROM, POS, ID, REF, ALT, QUAL, "PASS", INFO, "GT", gt_pred]
                    outline = '\t'.join([str(x) for x in outitems]) + '\n'
                    vcf_out.write(outline)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "The genotyping result has been written to the file:", vcf_out_file)
    print(current_time, "ALL DONE!")

if __name__ == "__main__":
        main()
