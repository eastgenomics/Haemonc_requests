import argparse
import pandas as pd
import pysam

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
        description="Information required to add GRCh38 lifted over variants"
    )
    parser.add_argument(
        "--csv",
        required=True,
        type=str,
        help="Path to CSV file to add GRCh38 lifted over variants to",
    )

    parser.add_argument(
        "--vcf",
        required=True,
        type=str,
        help="Path to VCF file which has been lifted over to GRCh38",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Name of output CSV file with GRCh38 lifted over variants added",
    )

    return parser.parse_args()


def read_vcf_to_df(vcf_file: str) -> pd.DataFrame:
    """
    Read in a VCF file into a dataframe.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the relevant columns from the VCF file.
    """
    vcf_in = pysam.VariantFile(vcf_file, "r")

    # Store required data from GRCh38 VCF in a list of dictionaries then
    # convert to a DataFrame
    records = []
    for rec in vcf_in:
        row = {
            "CHROM": str(rec.chrom),
            "POS": int(rec.pos),
            "REF": str(rec.ref),
            "ALT": ",".join(str(a) for a in rec.alts),
            "Genie_description": rec.info.get("Genie_description", None),
        }

        records.append(row)

    vcf_df = pd.DataFrame(records)

    return vcf_df


def merge_dataframes(
    b37_count_csv: pd.DataFrame, vcf_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge two dataframes on the common columns

    Parameters
    ----------
    b37_count_csv : pd.DataFrame
        DataFrame containing the b37 count CSV data.
    vcf_df : pd.DataFrame
        DataFrame containing the VCF data with GRCh38 lifted over.

    Returns
    -------
    merged_df : pd.DataFrame
        Merged DataFrame
    """
    merged_df = pd.merge(
        b37_count_csv,
        vcf_df,
        left_on="variant_description",
        right_on="Genie_description",
        how="left",
    )

    merged_df = merged_df.rename(
        columns={
            "CHROM": "chrom_grch38",
            "POS": "pos_grch38",
            "REF": "ref_grch38",
            "ALT": "alt_grch38",
            "chrom_vcf": "chrom_grch37",
            "pos_vcf": "pos_grch37",
            "ref_vcf": "ref_grch37",
            "alt_vcf": "alt_grch37",
        }
    )

    priority_cols = [
        "Hugo_Symbol",
        "chrom_grch38",
        "pos_grch38",
        "ref_grch38",
        "alt_grch38",
        "Genie_description",
        "HGVSc",
        "HGVSp",
        "RefSeq",
        "Consequence",
        "Variant_Classification",
    ]
    end_cols = [col for col in merged_df.columns if "_count" in col]
    middle_cols = [
        col for col in merged_df.columns if col not in priority_cols + end_cols
    ]
    merged_df = merged_df[priority_cols + middle_cols + end_cols]

    return merged_df


def main():
    args = parse_args()
    b37_count_csv = read_in_to_df(
        args.csv,
        sep=",",
        header=0,
    )
    b37_count_csv["variant_description"] = (
        b37_count_csv["Chromosome"].astype(str)
        + "_"
        + b37_count_csv["Start_Position"].astype(str)
        + "_"
        + b37_count_csv["Reference_Allele"].astype(str)
        + "_"
        + b37_count_csv["Tumor_Seq_Allele2"].astype(str)
    )

    vcf_df = read_vcf_to_df(args.vcf)
    merged_df = merge_dataframes(b37_count_csv, vcf_df)
    merged_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
