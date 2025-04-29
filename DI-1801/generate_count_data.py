import argparse
import pandas as pd


from utils import read_in_to_df, read_txt_file_to_list


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to query MAF file"
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Path to MAF file which includes patient and sample information",
    )

    parser.add_argument(
        "--haemonc_cancer_types",
        required=True,
        type=str,
        help=(
            "Path to file which lists haemonc cancer types we're interested in"
        ),
    )

    parser.add_argument(
        "--output", required=True, type=str, help="Name of output file"
    )

    return parser.parse_args()


def calculate_count_for_all_cancers(genie_data_with_sample_info):
    """
    Generate a new column with a count of the number of patients each variant
    has been seen in for all cancers, aggregating other relevant fields

    Parameters
    ----------
    genie_data_with_sample_info : pd.DataFrame
        Genie data to calculate counts from

    Returns
    -------
    overall_count_per_variant : pd.DataFrame
        dataframe with one row per variant, the count of that in patients
        and aggregated other fields
    """
    # Set data types so we don't end up with weird grouping due to formatting
    # differences
    genie_data_with_sample_info["Start_Position"] = (
        genie_data_with_sample_info["Start_Position"].astype("Int64")
    )
    genie_data_with_sample_info["Chromosome"] = (
        genie_data_with_sample_info["Chromosome"].astype(str).str.strip()
    )
    genie_data_with_sample_info["Reference_Allele"] = (
        genie_data_with_sample_info["Reference_Allele"].astype(str).str.strip()
    )
    genie_data_with_sample_info["Tumor_Seq_Allele2"] = (
        genie_data_with_sample_info["Tumor_Seq_Allele2"]
        .astype(str)
        .str.strip()
    )

    overall_count_per_variant = (
        genie_data_with_sample_info.drop_duplicates(
            subset=[
                "PATIENT_ID",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            keep="first",
        )
        .groupby(
            [
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ]
        )
        .agg(
            {
                "PATIENT_ID": "count",
                "Hugo_Symbol": lambda x: "&".join(sorted(x.dropna().unique())),
                "HGVSc": lambda x: "&".join(sorted(x.dropna().unique())),
                "HGVSp": lambda x: "&".join(sorted(x.dropna().unique())),
                "RefSeq": lambda x: "&".join(sorted(x.dropna().unique())),
                "Consequence": lambda x: "&".join(sorted(x.dropna().unique())),
                "Variant_Classification": lambda x: "&".join(
                    sorted(x.dropna().unique())
                ),
            }
        )
        .rename(columns={"PATIENT_ID": "overall_count"})
        .reset_index()
    )

    return overall_count_per_variant


def main():
    args = parse_args()
    genie_data_with_sample_info = read_in_to_df(
        args.input, header=0, dtype={"Chromosome": str}
    )

    overall_count_per_variant = calculate_count_for_all_cancers(
        genie_data_with_sample_info
    )

    overall_count_per_variant.to_csv(
        args.output,
        index=False,
    )

    # Overall count per cancer type
    haemonc_cancer_types = read_txt_file_to_list(args.haemonc_cancer_types)


if __name__ == "__main__":
    main()
