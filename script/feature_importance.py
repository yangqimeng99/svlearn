import argparse
import numpy as np
import polars as pl
import pandas as pd
from functools import reduce
from sklearn import utils
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from joblib import load


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

    condition1 = reduce(lambda a, b: a | b, [pl_df[col].is_null() for col in check_null])
    pl_df_filted = pl_df.filter(~condition1)

    condition2 = reduce(lambda a, b: a | b, [pl_df_filted[col].is_nan() for col in check_nan])
    pl_df_filted = pl_df_filted.filter(~condition2).cast(change_dtype_dict)

    return pl_df_filted, len(pl_df), len(pl_df_filted)


def main():
    parser = argparse.ArgumentParser(
        description="feature importance export tool: Support SVLearn 18/24 feature random forest model"
    )
    parser.add_argument("--model", required=True, help="model file")
    parser.add_argument("--train_set", required=True, help="training data set")
    parser.add_argument("--other_feature", action="store_true",
                        help="if this parameter is specified, the 24 feature scheme is used; otherwise, the 18 feature scheme is used")
    parser.add_argument("--out", required=True, help="feature importance result file name")
    parser.add_argument("--model_name", required=True, help="model name label")

    args = parser.parse_args()

    if args.other_feature:
        # 24 feature 情形
        categorical_features = [
            'sv_type', 'repeat_class', 'SDs_class', 'ref_GT', 'alt_GT', 'ref_FT', 'alt_FT'
        ]
        categories_per_feature = {
            'sv_type': [0, 1],
            'repeat_class': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            'SDs_class': [0, 1],
            'ref_GT': [0, 1, 2],
            'alt_GT': [0, 1, 2],
            'ref_FT': [0, 1, 2, 3, 4, 5, 6],
            'alt_FT': [0, 1, 2, 3, 4, 5, 6]
        }

        feature_dtype_change_dict = {
            "GT_true": pl.Int32,
            "sv_type": pl.Int32,
            "BP_lp_homref": pl.Float64,
            "BP_lp_het": pl.Float64,
            "BP_lp_homalt": pl.Float64,
            "ref_GT": pl.Int32,
            "alt_GT": pl.Int32,
            "PL": pl.Int32,
            "ref_FT": pl.Int32,
            "alt_FT": pl.Int32
        }
        check_null = ["BP_lp_homref", "BP_lp_het", "BP_lp_homalt"]
        check_nan = ["GT_true", "sv_type", "ref_GT", "alt_GT", "PL"]

    else:
        categorical_features = [
            'sv_type', 'repeat_class', 'SDs_class'
        ]
        categories_per_feature = {
            'sv_type': [0, 1],
            'repeat_class': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            'SDs_class': [0, 1]
        }
        feature_dtype_change_dict = {
            "GT_true": pl.Int32,
            "sv_type": pl.Int32,
            "BP_lp_homref": pl.Float64,
            "BP_lp_het": pl.Float64,
            "BP_lp_homalt": pl.Float64
        }
        check_null = ["BP_lp_homref", "BP_lp_het", "BP_lp_homalt"]
        check_nan = ["GT_true", "sv_type"]

    model = load(args.model)

    one_hot_encoders = {
        feature: OneHotEncoder(categories=[categories_per_feature[feature]], handle_unknown='ignore')
        for feature in categorical_features
    }
    transformers = [(feature, one_hot_encoders[feature], [feature]) for feature in categorical_features]
    column_transformer = ColumnTransformer(transformers=transformers, remainder='passthrough')

    train_set_df = pl.read_csv(args.train_set, separator='\t', has_header=True, null_values=["-"])
    train_set_filted_df, input_train_num, filtered_train_num = pre_data(train_set_df,
                                                                       feature_dtype_change_dict,
                                                                       check_null,
                                                                       check_nan)
    if "sv_id" in train_set_filted_df.columns:
        train_set_filted_df = train_set_filted_df.drop("sv_id")

    train_set_filted_df = utils.shuffle(train_set_filted_df)

    train_set_filted_df_pd = train_set_filted_df.to_pandas()

    data_y = train_set_filted_df_pd['GT_true'].values
    data_x = train_set_filted_df_pd.drop("GT_true", axis=1)

    column_transformer.fit_transform(data_x)
    transformed_feature_names = column_transformer.get_feature_names_out()

    feature_importances = model.feature_importances_

    original_feature_importances = {}
    for feature in categorical_features:
        original_feature_importances[feature] = 0.0

    for transformed_name, importance in zip(transformed_feature_names, feature_importances):
        found = False
        for original_feature in categorical_features:
            if transformed_name.startswith(original_feature):
                original_feature_importances[original_feature] += importance
                found = True
                break
        if not found:
            original_feature_importances[transformed_name] = importance

    result_rows = []
    for feat_name, imp in original_feature_importances.items():
        result_rows.append({
            "model_name": args.model_name,
            "feature": feat_name,
            "importance": imp
        })

    df_result = pd.DataFrame(result_rows)
    df_result.to_csv(args.out, index=False, sep='\t')
    print(f"[INFO] Feature importances saved to {args.out}")

if __name__ == "__main__":
    main()
