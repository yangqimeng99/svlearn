
import sys
import pysam
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

    precision_GT = len(genotype_df_pred_01_11_accuracy) / len(genotype_df_pred_01_11)
    benchmark_result['precision_GT'] = round(precision_GT, 4)

    genotype_df_true_01_11 = genotype_df.filter(pl.col('GT_true')!='./.').filter(pl.col('GT_true')!='0/0')
    genotype_df_true_01_11_accuracy = genotype_df_true_01_11.filter(pl.col('GT_true')==pl.col('GT_pred'))
    recall_GT = len(genotype_df_true_01_11_accuracy) / len(genotype_df_true_01_11)
    benchmark_result['recall_GT'] = round(recall_GT, 4)

    f1_GT = ( 2* precision_GT * recall_GT) / (precision_GT + recall_GT)
    benchmark_result['f1_GT'] = round(f1_GT, 4)

    genotype_df_pred_01_11_accuracy_2 = (genotype_df_pred_01_11
                                         .filter(pl.col('GT_true')!='0/0')
                                         .filter(pl.col('GT_true')!='./.'))
    precision = len(genotype_df_pred_01_11_accuracy_2) / len(genotype_df_pred_01_11)
    benchmark_result['precision'] = round(precision, 4)

    genotype_df_true_01_11_accuracy_2 = (genotype_df_true_01_11
                                         .filter(pl.col('GT_pred')!='0/0')
                                         .filter(pl.col('GT_pred')!='./.'))
    recall = len(genotype_df_true_01_11_accuracy_2) / len(genotype_df_true_01_11)
    benchmark_result['recall'] = round(recall, 4)

    f1 = ( 2* precision * recall) / (precision + recall)
    benchmark_result['f1'] = round(f1, 4)
    
    genotype_df_00 = genotype_df.filter(pl.col('GT_true')=='0/0').filter(pl.col('GT_pred')!='./.')
    genotype_df_00_accuracy = genotype_df_00.filter(pl.col('GT_true')==pl.col('GT_pred'))
    conc_00 = len(genotype_df_00_accuracy) / len(genotype_df_00)
    benchmark_result['conc_00'] = round(conc_00, 4)

    genotype_df_01 = genotype_df.filter(pl.col('GT_true')=='0/1').filter(pl.col('GT_pred')!='./.')
    genotype_df_01_accuracy = genotype_df_01.filter(pl.col('GT_true')==pl.col('GT_pred'))
    conc_01 = len(genotype_df_01_accuracy) / len(genotype_df_01)
    benchmark_result['conc_01'] = round(conc_01, 4)

    genotype_df_11 = genotype_df.filter(pl.col('GT_true')=='1/1').filter(pl.col('GT_pred')!='./.')
    genotype_df_11_accuracy = genotype_df_11.filter(pl.col('GT_true')==pl.col('GT_pred'))
    conc_11 = len(genotype_df_11_accuracy) / len(genotype_df_11)
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
    parser = argparse.ArgumentParser(description="Benchmark the genotyped SV set based SV ID.")
    parser.add_argument("-b", "--base_set", type=str, required=True, help="True sv vcf set", metavar="file")
    parser.add_argument("-c", "--call_set", type=str, required=True, help="Predict sv vcf set", metavar="file")
    parser.add_argument("-n1", "--base_sample_name", type=str, required=True, help="Sample name in true sv vcf set", metavar="str")
    parser.add_argument("-n2", "--call_sample_name", type=str, required=True, help="Sample name in predict sv vcf set", metavar="str")
    parser.add_argument("-o", "--out", type=str, default="benchmark.result.tsv", help="The out of bechmark result, Default: benchmark.result.tsv", metavar="file")

    parsed_args = parser.parse_args(args=args)
    base_set_file = parsed_args.base_set
    call_set_file = parsed_args.call_set
    base_sample_name = parsed_args.base_sample_name
    call_sample_name = parsed_args.call_sample_name
    benchmark_out_file = parsed_args.out

    base_set = pysam.VariantFile(base_set_file, "r")
    call_set = pysam.VariantFile(call_set_file, "r")
    base_num_records = len(list(base_set.fetch()))
    call_num_records = len(list(call_set.fetch()))
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "The True sv number is:", base_num_records)
    print(current_time, "The Predict sv number is:", call_num_records)
    base_set.close()
    call_set.close()

    base_set = pysam.VariantFile(base_set_file, "r")
    call_set = pysam.VariantFile(call_set_file, "r")
    base_call_set = {}

    gt_dict = {
        (0, 1): "0/1",
        (1, 0): "0/1",
        (0, 0): "0/0",
        (1, 1): "1/1",
        (None,): "./.",
    }

    for call_record in call_set:    
        call_gt = gt_dict.get(call_record.samples[call_sample_name]['GT'], "./.")
        base_call_set[call_record.id] = [call_gt]

    for base_record in base_set:
        if base_record.id in base_call_set:
            true_gt = gt_dict.get(base_record.samples[base_sample_name]['GT'], "./.")
            base_call_set[base_record.id].insert(0, true_gt)
        else:
            base_call_set[base_record.id] = [true_gt, "./."]

    base_set.close()
    call_set.close()

    sv_id = list(base_call_set.keys())
    true_GT = [value[0] for value in base_call_set.values()]
    pred_GT = [value[1] for value in base_call_set.values()]

    geno_df = pl.DataFrame({
        'sv_id': sv_id,
        'GT_true': true_GT,
        'GT_pred': pred_GT
    })

    geno_df = geno_df.filter(pl.col('GT_true')!='./.')

    if len(geno_df) == 0:
        print("No sv can be genotyped.")
        sys.exit(1)
    elif len(geno_df) < base_num_records:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(current_time, "The filtered true sv number is:", len(geno_df))
        print(current_time, "The discarded true sv number is:", base_num_records - len(geno_df))

    get_benchmark_info(geno_df, benchmark_out_file)

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(current_time, "Benchmark result has been saved in", benchmark_out_file)
    print(current_time, "Benchmark finished.")

if __name__ == "__main__":
    main()


