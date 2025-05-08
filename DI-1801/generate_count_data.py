import argparse
import pandas as pd


from utils import read_in_to_df, read_txt_file_to_list

VARIANT_COLUMNS = [
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
]


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to generate counts from MAF file"
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
    Generate a dataframe with one row per unique variant, a count of how
    many patients it occurs in (all cancers) with other relevant fields
    aggregated

    Parameters
    ----------
    genie_data_with_sample_info : pd.DataFrame
        Genie data with correct datatypes set to calculate counts from

    Returns
    -------
    overall_count_per_variant_reordered : pd.DataFrame
        dataframe with one row per variant, the count of that in patients
        and aggregated other fields
    """
    # Drop duplicates per variant, group by each variant and count how many
    # patients have that variant, aggregating other fields with an & char
    subset_fields = ["PATIENT_ID"] + VARIANT_COLUMNS
    overall_count_per_variant = (
        genie_data_with_sample_info.drop_duplicates(
            subset=subset_fields,
            keep="first",
        )
        .groupby(VARIANT_COLUMNS)
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
        .rename(columns={"PATIENT_ID": "all_cancers_count"})
        .reset_index()
    )

    # Reorder columns so we always keep these first
    first_cols = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "HGVSc",
        "HGVSp",
        "RefSeq",
        "Consequence",
        "Variant_Classification",
    ]
    overall_count_per_variant_reordered = overall_count_per_variant[
        first_cols
        + [
            col
            for col in overall_count_per_variant.columns
            if col not in first_cols
        ]
    ]

    return overall_count_per_variant_reordered


def calculate_count_for_all_haemonc_cancers(
    genie_data_with_sample_info, haemonc_cancer_types
):
    """
    Create df with one row per unique variant and a count of each variant in
    haemonc cancers

    Parameters
    ----------
    genie_data_with_sample_info : pd.DataFrame
        Genie data to calculate counts from
    haemonc_cancer_types : list
        List of haemonc cancer types (strings) to subset to (must match the
        names of the cancer types in the Genie file)

    Returns
    -------
    haemonc_cancer_counts : pd.DataFrame
        Count data for each variant in haemonc cancers
    """
    # Create multiindex of all unique variants in the dataset
    all_variant_index = genie_data_with_sample_info.drop_duplicates(
        subset=VARIANT_COLUMNS
    )[VARIANT_COLUMNS]

    # Create index from the four variant columns
    all_variant_index = all_variant_index.set_index(
        VARIANT_COLUMNS
    ).sort_index()

    # Subset to just haemonc cancer types
    haemonc_cancer_type_rows = genie_data_with_sample_info[
        genie_data_with_sample_info["CANCER_TYPE"].isin(haemonc_cancer_types)
    ]

    subset_fields = ["PATIENT_ID"] + VARIANT_COLUMNS
    # Group by patient and variant and count
    haemonc_cancer_counts = (
        haemonc_cancer_type_rows.drop_duplicates(subset=subset_fields)
        .groupby(VARIANT_COLUMNS)
        .agg(haemonc_cancers_count=("PATIENT_ID", "count"))
    )

    # Use the index made earlier to set any variants which aren't present in
    # any haemonc cancers to a count of zero
    haemonc_cancer_counts = haemonc_cancer_counts.reindex(
        all_variant_index.index, fill_value=0
    ).reset_index()

    return haemonc_cancer_counts


def calculate_count_per_haemonc_cancer_type(
    genie_data_with_sample_info, haemonc_cancer_types
):
    """
    Generate df with one row per unique variant and a count of each variant
    in each specific haemonc cancer type

    Parameters
    ----------
    genie_data_with_sample_info : pd.DataFrame
        Genie data to calculate counts from
    haemonc_cancer_types : list
        List of haemonc cancer types (strings) to subset to (must match the
        names of the cancer types in the Genie file)

    Returns
    -------
    per_cancer_haemonc_subset : pd.DataFrame
        dataframe with one row per variant and columns with counts of that
        variant for each haemonc cancer type
    """
    subset_fields = ["PATIENT_ID"] + VARIANT_COLUMNS + ["CANCER_TYPE"]
    groupby_fields = VARIANT_COLUMNS + ["CANCER_TYPE"]
    per_cancer_counts = (
        genie_data_with_sample_info.drop_duplicates(
            subset=subset_fields,
            keep="first",
        )
        .groupby(groupby_fields)
        .size()
        .reset_index(name="variant_count")
    )

    pivoted_per_cancer_counts = per_cancer_counts.pivot_table(
        index=VARIANT_COLUMNS,
        columns="CANCER_TYPE",
        values="variant_count",
        fill_value=0,
    ).reset_index()

    # Remove any columns for counts which aren't in the haemonc cancer types
    missing_cancer_types = [
        col
        for col in haemonc_cancer_types
        if col not in pivoted_per_cancer_counts.columns
    ]
    if missing_cancer_types:
        print(
            "Warning: the following haemonc cancer types are missing from"
            f" the data: {', '.join(missing_cancer_types)}"
        )
    existing_cancer_types = [
        col for col in haemonc_cancer_types if col in pivoted_per_cancer_counts
    ]
    all_cols_to_keep = VARIANT_COLUMNS + existing_cancer_types

    # Subset to just counts for haemonc cancer types
    per_cancer_haemonc_subset = pivoted_per_cancer_counts[all_cols_to_keep]
    # Add '_count' suffix to each count column
    per_cancer_haemonc_subset.columns = [
        (col if col in VARIANT_COLUMNS else f"{col}_count")
        for col in per_cancer_haemonc_subset.columns
    ]

    return per_cancer_haemonc_subset


def merge_counts_together(df1, df2):
    """
    Merge the counts together to get a single dataframe with multiple counts
    (all have one row per unique variant so will have same number of rows)

    Parameters
    ----------
    df1 : pd.DataFrame
        First dataframe to merge
    df2 : pd.DataFrame
        Second dataframe to merge

    Returns
    -------
    merged_counts : pd.DataFrame
        Merged dataframe with all counts
    """
    merged_counts = pd.merge(
        df1,
        df2,
        on=VARIANT_COLUMNS,
        how="left",
    )

    return merged_counts


def main():
    args = parse_args()
    genie_data_with_sample_info = read_in_to_df(
        args.input,
        header=0,
        dtype={
            "Entrez_Gene_Id": "Int64",
            "Start_Position": "Int64",
        },
        converters={
            col: lambda x: x.strip() if isinstance(x, str) else x
            for col in ["Chromosome", "Reference_Allele", "Tumor_Seq_Allele2"]
        },
    )
    all_cancers_count = calculate_count_for_all_cancers(
        genie_data_with_sample_info
    )

    haemonc_cancer_types = read_txt_file_to_list(args.haemonc_cancer_types)
    all_haemonc_cancers_count = calculate_count_for_all_haemonc_cancers(
        genie_data_with_sample_info, haemonc_cancer_types
    )
    count_per_haemonc_cancer_type = calculate_count_per_haemonc_cancer_type(
        genie_data_with_sample_info, haemonc_cancer_types
    )

    all_cancers_and_haemonc_cancers = merge_counts_together(
        all_cancers_count, all_haemonc_cancers_count
    )
    merged_all = merge_counts_together(
        all_cancers_and_haemonc_cancers, count_per_haemonc_cancer_type
    )
    # Write to CSV
    merged_all.to_csv(
        args.output,
        index=False,
    )


if __name__ == "__main__":
    main()
