import sys
import argparse
import polars as pl
from datetime import datetime

def get_benchmark_info(genotype_df, benchmark_out_file):
    benchmark_result = {}
    sv_set_num = len(genotype_df)
    benchmark_result['sv_set_number'] = sv_set_num

    genotyped_df = genotype_df.filter(pl.col("GT_pred") != "./.")
    genotyped_sv_num = len(genotyped_df)
    benchmark_result['genotyped_sv_number'] = genotyped_sv_num

    genotype_rate = (genotyped_sv_num / sv_set_num)
    benchmark_result['genotype_rate'] = round(genotype_rate, 4)

    genotype_df_right = genotype_df.filter(pl.col('GT_pred')==pl.col('GT_true'))
    accuracy_genotyped_sv_number = len(genotype_df_right)
    benchmark_result['accuracy_genotyped_sv_number'] = accuracy_genotyped_sv_number

    genotype_df_pred_01_11 = genotype_df.filter(pl.col('GT_pred')!='./.').filter(pl.col('GT_pred')!='0/0')
    genotype_df_pred_01_11_accuracy = genotype_df_pred_01_11.filter(pl.col('GT_true')==pl.col('GT_pred'))

    precison_GT = len(genotype_df_pred_01_11_accuracy) / len(genotype_df_pred_01_11) if len(genotype_df_pred_01_11) > 0 else 0.0
    benchmark_result['precison_GT'] = round(precison_GT, 4)

    genotype_df_true_01_11 = genotype_df.filter(pl.col('GT_true')!='./.').filter(pl.col('GT_true')!='0/0')
    genotype_df_true_01_11_accuracy = genotype_df_true_01_11.filter(pl.col('GT_true')==pl.col('GT_pred'))
    recall_GT = len(genotype_df_true_01_11_accuracy) / len(genotype_df_true_01_11) if len(genotype_df_true_01_11) > 0 else 0.0
    benchmark_result['recall_GT'] = round(recall_GT, 4)

    if precison_GT + recall_GT > 0:
        f1_GT = ( 2* precison_GT * recall_GT) / (precison_GT + recall_GT)
    else:
        f1_GT = 0.0
    benchmark_result['f1_GT'] = round(f1_GT, 4)

    genotype_df_pred_01_11_accuracy_2 = (genotype_df_pred_01_11
                                         .filter(pl.col('GT_true')!='0/0')
                                         .filter(pl.col('GT_true')!='./.'))
    precison = len(genotype_df_pred_01_11_accuracy_2) / len(genotype_df_pred_01_11) if len(genotype_df_pred_01_11) > 0 else 0.0
    benchmark_result['precison'] = round(precison, 4)

    genotype_df_true_01_11_accuracy_2 = (genotype_df_true_01_11
                                         .filter(pl.col('GT_pred')!='0/0')
                                         .filter(pl.col('GT_pred')!='./.'))
    recall = len(genotype_df_true_01_11_accuracy_2) / len(genotype_df_true_01_11) if len(genotype_df_true_01_11) > 0 else 0.0
    benchmark_result['recall'] = round(recall, 4)

    if precison + recall > 0:
        f1 = ( 2* precison * recall) / (precison + recall)
    else:
        f1 = 0.0
    benchmark_result['f1'] = round(f1, 4)

    genotype_df_00 = genotype_df.filter(pl.col('GT_true')=='0/0').filter(pl.col('GT_pred')!='./.')
    if len(genotype_df_00) > 0:
        genotype_df_00_accuracy = genotype_df_00.filter(pl.col('GT_true')==pl.col('GT_pred'))
        conc_00 = len(genotype_df_00_accuracy) / len(genotype_df_00)
    else:
        conc_00 = 0.0
    benchmark_result['conc_00'] = round(conc_00, 4)

    genotype_df_01 = genotype_df.filter(pl.col('GT_true')=='0/1').filter(pl.col('GT_pred')!='./.')
    if len(genotype_df_01) > 0:
        genotype_df_01_accuracy = genotype_df_01.filter(pl.col('GT_true')==pl.col('GT_pred'))
        conc_01 = len(genotype_df_01_accuracy) / len(genotype_df_01)
    else:
        conc_01 = 0.0
    benchmark_result['conc_01'] = round(conc_01, 4)

    genotype_df_11 = genotype_df.filter(pl.col('GT_true')=='1/1').filter(pl.col('GT_pred')!='./.')
    if len(genotype_df_11) > 0:
        genotype_df_11_accuracy = genotype_df_11.filter(pl.col('GT_true')==pl.col('GT_pred'))
        conc_11 = len(genotype_df_11_accuracy) / len(genotype_df_11)
    else:
        conc_11 = 0.0
    benchmark_result['conc_11'] = round(conc_11, 4)

    wgc = (conc_00 + conc_01 + conc_11) / 3
    benchmark_result['wgc'] = round(wgc, 4)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "Verify model effect:")
    with open(benchmark_out_file, "w") as file:
        for key, value in benchmark_result.items():
            print(f"{key}\t{value}")
            file.write(f"{key}  {value}\n")

def main(args=None):
    parser = argparse.ArgumentParser(description="Benchmark the genotyped SV set based on a simple TSV input.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input TSV file containing sv_id, GT_true, GT_pred")
    parser.add_argument("-o", "--out", type=str, default="benchmark.result.tsv", help="Output benchmark result file")
    parsed_args = parser.parse_args(args=args)

    input_file = parsed_args.input
    benchmark_out_file = parsed_args.out

  
    geno_df = pl.read_csv(input_file, separator='\t')

    geno_df = geno_df.filter(pl.col('GT_true')!='./.')

    if len(geno_df) == 0:
        print("No sv can be genotyped.")
        sys.exit(1)

    get_benchmark_info(geno_df, benchmark_out_file)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "Benchmark result has been saved in", benchmark_out_file)
    print(current_time, "Benchmark finished.")

if __name__ == "__main__":
    main()
