# %%
import os
import sys
import argparse
import subprocess
import numpy as np
import polars as pl
import pandas as pd

from functools import reduce
from datetime import datetime

from sklearn import utils
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV
from joblib import dump


def get_output_file_path(output_dir_path, file_name):
	return os.path.join(output_dir_path, file_name)


def pre_data(pl_df, change_dtype_dict, check_null, check_nan):
	pl_df = pl_df.with_columns(
		pl.when(pl_df['GT_true'] == "0/0").then(0)
		.when(pl_df['GT_true'] == "0/1").then(1)
		.when(pl_df['GT_true'] == "1/1").then(2)
		.otherwise(pl.lit(np.nan)).alias("GT_true")
		)
	
	pl_df = pl_df.with_columns(
		pl.when(pl_df['sv_type'] == "INS").then(0)
		.when(pl_df['sv_type'] == "DEL").then(1)
		.otherwise(pl.lit(np.nan)).alias("sv_type")
		)

	condition1 = reduce(lambda a, b: a | b, [pl_df[col] == "-" for col in check_null])
	pl_df_filted = pl_df.filter(~condition1)

	condition2 = reduce(lambda a, b: a | b, [pl_df_filted[col].is_nan() for col in check_nan])
	pl_df_filted = pl_df_filted.filter(~condition2).cast(change_dtype_dict)

	return pl_df_filted, len(pl_df), len(pl_df_filted)



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

	precison_GT = len(genotype_df_pred_01_11_accuracy) / len(genotype_df_pred_01_11)
	benchmark_result['precison_GT'] = round(precison_GT, 4)

	genotype_df_true_01_11 = genotype_df.filter(pl.col('GT_true')!='./.').filter(pl.col('GT_true')!='0/0')
	genotype_df_true_01_11_accuracy = genotype_df_true_01_11.filter(pl.col('GT_true')==pl.col('GT_pred'))
	recall_GT = len(genotype_df_true_01_11_accuracy) / len(genotype_df_true_01_11)
	benchmark_result['recall_GT'] = round(recall_GT, 4)

	f1_GT = ( 2* precison_GT * recall_GT) / (precison_GT + recall_GT)
	benchmark_result['f1_GT'] = round(f1_GT, 4)

	genotype_df_pred_01_11_accuracy_2 = (genotype_df_pred_01_11
										 .filter(pl.col('GT_true')!='0/0')
										 .filter(pl.col('GT_true')!='./.'))
	precison = len(genotype_df_pred_01_11_accuracy_2) / len(genotype_df_pred_01_11)
	benchmark_result['precison'] = round(precison, 4)

	genotype_df_true_01_11_accuracy_2 = (genotype_df_true_01_11
										 .filter(pl.col('GT_pred')!='0/0')
										 .filter(pl.col('GT_pred')!='./.'))
	recall = len(genotype_df_true_01_11_accuracy_2) / len(genotype_df_true_01_11)
	benchmark_result['recall'] = round(recall, 4)

	f1 = ( 2* precison * recall) / (precison + recall)
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

	print("Verify model effect:")
	with open(benchmark_out_file, "w") as file:
		for key, value in benchmark_result.items():
			print(f"{key}\t{value}")
			file.write(f"{key}  {value}\n")


# %%
def main(args=None):
	parser = argparse.ArgumentParser(description="Run paragraph in ref and alt bam.")
	parser.add_argument("--train_set", type=str, required=True, help="Training data set", metavar="file")
	parser.add_argument("--val_set", type=str, default=None, help="Verification data set", metavar="file")
	parser.add_argument('--other_feature', type=str, default=None, help='optional:{paragraph, other_feature}', metavar="str")
	parser.add_argument('--train_model', type=str, required=True, help='optional: {RandomForest}', metavar="str")
	parser.add_argument('-t', '--threads', type=str, default="1", help='The number of threads to use', metavar="int")
	parser.add_argument('-o', '--out', type=str, default="trainingModel_out", help="The output dir name, Default: trainingModel_out", metavar="dir")

	parsed_args = parser.parse_args(args=args)
	train_set = parsed_args.train_set
	val_set = parsed_args.val_set
	other_feature = parsed_args.other_feature
	train_model = parsed_args.train_model
	threads = parsed_args.threads
	output_dir = parsed_args.out

	other_feature = None if other_feature == "None" else other_feature
	val_set = None if val_set == "None" else val_set
	
	if other_feature is None:
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "No other feature input")
	elif other_feature not in ["paragraph", "other_feature"]:
		print("Error: please input the correct feature type!")
		exit(1)

	threads = int(threads)
	output_dir_path = os.path.abspath(output_dir)
	os.makedirs(output_dir_path, exist_ok=True)

	train_set_df = pl.read_csv(train_set, separator='\t', has_header=True)

	if other_feature == "paragraph":
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "paragraph feature input")
		feature_dtype_change_dict = {
			"GT_true":pl.Int32,
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
		check_nan = ["GT_true", "sv_type", "ref_GT", "alt_GT", "PL"]
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
	elif other_feature is None:
		feature_dtype_change_dict = {
			"GT_true":pl.Int32,
			"sv_type":pl.Int32,
			"BP_lp_homref":pl.Float64, 
			"BP_lp_het":pl.Float64, 
			"BP_lp_homalt":pl.Float64
			}
		check_null = ["BP_lp_homref", "BP_lp_het", "BP_lp_homalt"]
		check_nan = ["GT_true", "sv_type"]
		categorical_features = ['sv_type', 'repeat_class', 'SDs_class']
		categories_per_feature = {
			'sv_type': [0, 1],
			'repeat_class': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
			'SDs_class': [0, 1]
		}
	else:
		print("Other features are not currently supported")
		exit(1)

	train_set_filted_df, input_train_num, filtered_train_num, = pre_data(train_set_df, feature_dtype_change_dict, check_null, check_nan)
	train_set_filted_df = train_set_filted_df.drop("sv_id")
	train_set_filted_df = utils.shuffle(train_set_filted_df)

	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "The number of training data entries entered is:", input_train_num)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "The number of training data entries after filtering is:", filtered_train_num)
	print(current_time, "The number of Feature entries entered is:", train_set_filted_df.shape[1]-1)

	train_set_filted_df_pd = train_set_filted_df.to_pandas()
	data_y = train_set_filted_df_pd['GT_true'].values
	data_x = train_set_filted_df_pd.drop("GT_true", axis=1)

	one_hot_encoders = {
		feature: OneHotEncoder(categories=[categories_per_feature[feature]], handle_unknown='ignore')
		for feature in categorical_features
	}
	transformers = [(feature, one_hot_encoders[feature], [feature]) for feature in categorical_features]
	column_transformer = ColumnTransformer(transformers=transformers, remainder='passthrough')

	train_X_encoded = column_transformer.fit_transform(data_x)

	skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
	if train_model == "RandomForest":

		param_grid = {
			'n_estimators': [100, 200, 300, 400, 500, 750, 1000],
			'max_depth': [None, 10, 20, 30, 50, 100, 200, 300],
			'min_samples_split': [2, 5, 10],
			'min_samples_leaf': [1, 2, 4],
			'bootstrap': [True, False],
			'class_weight': [None, 'balanced']
		}
		
		clf = HalvingGridSearchCV(RandomForestClassifier(), param_grid, cv=skf, n_jobs=threads, scoring='accuracy', verbose=2)
		clf.fit(train_X_encoded, data_y)

		best_params = clf.best_params_
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, f"Best parameters: {best_params}")
		results_train_df = pd.DataFrame(clf.cv_results_)
		RandomForest_grid_search_results_out = get_output_file_path(output_dir_path, "RandomForest_grid_search_results.tsv")
		results_train_df.to_csv(RandomForest_grid_search_results_out, sep='\t', index=False)

		model = RandomForestClassifier(**best_params, n_jobs=threads)
		model.fit(train_X_encoded, data_y)
		RandomForest_model_out = get_output_file_path(output_dir_path, "RandomForest_model.joblib")
		dump(model, RandomForest_model_out, compress=5)

	else:
		print("Please select a supported model for training: RandomForest")
		sys.exit(1)

	if val_set is None:
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "ALL DONE!")
		sys.exit(0)
	else:
		val_set_df = pl.read_csv(val_set, separator='\t', has_header=True)
		val_set_filted_id_df, input_val_num, filtered_val_num, = pre_data(val_set_df, feature_dtype_change_dict, check_null, check_nan)
		val_set_filted_df = val_set_filted_id_df.drop("sv_id")
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "The number of validation data entries entered is:", input_val_num)
		print(current_time, "The number of validation data entries after filtering is:", filtered_val_num)

		val_set_filted_df_pd = val_set_filted_df.to_pandas()
		val_y = val_set_filted_df_pd['GT_true'].values
		val_x = val_set_filted_df_pd.drop("GT_true", axis=1)
		val_x_encoded = column_transformer.fit_transform(val_x)
		pred_y_val = model.predict(val_x_encoded)
		accuracy = accuracy_score(val_y, pred_y_val)

		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "The validation data accuracy is:", accuracy)

		val_set_filted_id_df = (val_set_filted_id_df
						.with_columns(pl.Series("GT_pred", pred_y_val))
						.select(["sv_id", "GT_pred"]))
		
		val_out_df = (val_set_df
					  .join(val_set_filted_id_df, on='sv_id', how='left', validate='1:1')
					  .select(["sv_id", "GT_true", "GT_pred"])
					  .with_columns(
						  pl.when(pl.col('GT_pred') == 0).then(pl.lit("0/0"))
						  .when(pl.col('GT_pred') == 1).then(pl.lit("0/1"))
						  .when(pl.col('GT_pred') == 2).then(pl.lit("1/1"))
						  .otherwise(pl.lit("./.")).alias("GT_pred")))

		val_out_df.write_csv(get_output_file_path(output_dir_path, "val_set_pred.tsv"), separator="\t")
		val_result_out = get_output_file_path(output_dir_path, "val_result.txt")
		get_benchmark_info(val_out_df, val_result_out)
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "ALL DONE!")

if __name__ == "__main__":
	main()
