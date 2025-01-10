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
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

# %%
def get_output_file_path(output_dir_path, file_name):
	return os.path.join(output_dir_path, file_name)

# %%
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
	
	# condition1 = reduce(lambda a, b: a | b, [pl_df[col] == "-" for col in check_null])
	condition1 = reduce(lambda a, b: a | b, [pl_df[col].is_null() for col in check_null])
	pl_df_filted = pl_df.filter(~condition1)

	condition2 = reduce(lambda a, b: a | b, [pl_df_filted[col].is_nan() for col in check_nan])
	pl_df_filted = pl_df_filted.filter(~condition2).cast(change_dtype_dict)

	return pl_df_filted, len(pl_df), len(pl_df_filted)

# %%
def main(args=None):
	parser = argparse.ArgumentParser(description="Choose Machine learning model.")
	parser.add_argument("--train_set", type=str, required=True, help="Training data set")
	parser.add_argument('--other_feature', type=str, default=None, help='{paragraph, other_feature}')
	parser.add_argument('--threads', type=str, default="1", help='The number of threads to use.')
	parser.add_argument("--out", "-o", type=str, default="trainingModel_out", help="The output dir name.")

	parsed_args = parser.parse_args(args=args)
	train_set = parsed_args.train_set
	other_feature = parsed_args.other_feature
	threads = parsed_args.threads
	output_dir = parsed_args.out

	other_feature = None if other_feature == "None" else other_feature
	
	if other_feature is None:
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, "No other feature input")
	elif other_feature not in ["paragraph", "other_feature"]:
		print("Error: please input the correct feature type!")
		exit(1)

	threads = int(threads)
	output_dir_path = os.path.abspath(output_dir)
	os.makedirs(output_dir_path, exist_ok=True)

	train_set_df = pl.read_csv(train_set, separator='\t', has_header=True, null_values=["-"])

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
	models = {
		"Random Forest": RandomForestClassifier(random_state=42, n_jobs=threads),
		"K-Nearest Neighbors": KNeighborsClassifier(n_jobs=threads),
		"Naive Bayes": GaussianNB(),
		"Logistic Regression": LogisticRegression(random_state=42, n_jobs=threads),
		"Gradient Boosting": GradientBoostingClassifier(random_state=42),
		"SVM": SVC(random_state=42)
	}

	results_all = {}
	results_mean = {}
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, "Start choose model...")
	# Perform 10-fold cross-validation and calculate average accuracy
	for name, model in models.items():
		pipeline = make_pipeline(StandardScaler(), model)
		cv_scores = cross_val_score(pipeline, train_X_encoded, data_y, cv=skf, scoring='accuracy')
		results_all[name] = cv_scores
		results_mean[name] = np.mean(cv_scores)
		current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print(current_time, f"{name}: Average 10-Fold Accuracy = {np.mean(cv_scores):.4f}")

	# Select the model with the highest accuracy
	best_model = max(results_mean, key=results_mean.get)
	current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	print(current_time, f"Best Model: {best_model}")
	
	choose_model_train_out = get_output_file_path(output_dir_path, "choose_model_train_cv_scores.tsv")
	with open(choose_model_train_out, 'w') as f:
		f.write('\t'.join(results_all.keys()) + '\n')
		for row in zip(*results_all.values()):
			f.write('\t'.join(map(str, row)) + '\n')

if __name__ == "__main__":
	main()
