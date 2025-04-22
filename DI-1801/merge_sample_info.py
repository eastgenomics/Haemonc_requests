import argparse
import pandas as pd

from utils import read_in_to_df


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to merge MAF file with sample info"
    )

    parser.add_argument(
        "--input_maf", required=True, type=str, help="Path to MAF file"
    )

    parser.add_argument(
        "--clinical_info",
        required=True,
        type=str,
        help=(
            "Path to file which lists sample IDs, patient IDs and cancer types"
        ),
    )

    parser.add_argument(
        "--output", required=True, type=str, help="Name of output merged file"
    )

    return parser.parse_args()


def merge_dataframes(
    df1, df2, left_on: str, right_on: str, how: str
) -> pd.DataFrame:
    """
    Merge two dataframes on common columns

    Parameters
    ----------
    df1 : pd.DataFrame
        first dataframe to merge
    df2 : pd.DataFrame
        second dataframe to merge
    left_on : str
        name of column which is common to both dataframes in first df
    right_on : str
        name of column which is common to both dataframes in second df

    Returns
    -------
    pd.DataFrame
        merged dataframe
    """
    try:
        merged_df = pd.merge(
            df1, df2, left_on=left_on, right_on=right_on, how=how
        )
        return merged_df
    except Exception as err:
        print(f"Error merging dataframes: {err}")
        return None


def perform_merge_checks(genie_data, genie_with_sample_ids):
    """
    Perform checks on merged dataframe

    Parameters
    ----------
    genie_data: pd.DataFrame
        original dataframe
    genie_with_sample_ids : pd.DataFrame
        merged dataframe

    Raises
    ------
    ValueError
        If the number of rows in the input and output data do not match
    """
    print(f"\nRows in input data: {len(genie_data)}")

    print(f"Rows in output data: {len(genie_with_sample_ids)}")
    if len(genie_data) != len(genie_with_sample_ids):
        raise ValueError(
            "Number of rows in input and output data do not match"
        )
    print("\nRows found with no patient ID after merging:")
    rows_with_na = genie_with_sample_ids[
        genie_with_sample_ids["PATIENT_ID"].isna()
    ]
    if not rows_with_na.empty:
        print(rows_with_na)
    else:
        print("None")


def main():
    args = parse_args()
    genie_data = read_in_to_df(
        args.input_maf,
        header=0,
    )
    sample_cancer_info = read_in_to_df(
        args.clinical_info,
        skiprows=4,
        header=0,
    )
    genie_with_sample_ids = merge_dataframes(
        genie_data,
        sample_cancer_info,
        left_on="Tumor_Sample_Barcode",
        right_on="SAMPLE_ID",
        how="inner",
    )
    perform_merge_checks(genie_data, genie_with_sample_ids)
    genie_with_sample_ids.to_csv(
        args.output,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
