import argparse
import pandas as pd
import os
import subprocess

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
        description="Information required to add GRCh38 coords"
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help=(
            "Path to Genie file with one row per variant and aggregated"
            " information"
        ),
    )

    parser.add_argument(
        "--chain",
        required=True,
        type=str,
        help=(
            "Path to liftover UCSC chain file (e.g. hg19ToHg38.over.chain.gz)"
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Path to output file with GRCh38 coordinates",
    )

    return parser.parse_args()


def split_liftover_bed_columns(df, colname):
    """
    Split the coordinates in the format "chr:start-end" into separate columns

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe with column to split
    colname : str
        name of the column we're splitting

    Returns
    -------
    df : pd.DataFrame
        pandas dataframe with columns split
    """
    # Split by colon ":" to separate chromosome and start-end
    parts = df[colname].str.split(":", expand=True)
    # Remove chr prefix
    df["GRCh38_chr"] = parts[0].str.removeprefix("chr")

    # Split start-end to separate columns
    start_end = parts[1].str.split("-", expand=True)
    df["GRCh38_start"] = start_end[0].astype(int)
    df["GRCh38_end"] = start_end[1].astype(int)

    df.drop(
        columns=["GRCh38_position"],
        inplace=True,
    )

    return df


def create_positions_file_from_count_df(genie_count_df):
    """
    Create and write out a positions file for liftover

    Parameters
    ----------
    genie_count_df : pd.DataFrame
        dataframe with one row per variant and aggregate counts
    """
    bed_df = genie_count_df[["Chromosome", "Start_Position", "End_Position"]]

    # Add 'chr' prefix which is needed for liftover
    bed_df["Chromosome"] = (
        bed_df["Chromosome"]
        .astype(str)
        .apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
    )

    # Format as chrX:start-end
    bed_df["positions"] = bed_df.apply(
        lambda row: f"{row['Chromosome']}:{row['Start_Position']}-{row['End_Position']}",
        axis=1,
    )

    bed_df[["positions"]].to_csv(
        "liftover_positions.txt",
        sep="\t",
        header=False,
        index=False,
    )


def run_liftover(chain_file):
    """
    Runs the UCSC liftOver tool on a BED file using the specified chain file.
    Assumes 'liftover_positions.txt' exists in the current directory.

    Parameters
    ----------
    chain_file : str
        Path to the liftover chain file (e.g. 'hg19ToHg38.over.chain.gz')
    """
    if not os.path.isfile("./liftOver"):
        raise FileNotFoundError(
            "liftOver binary not found at './liftOver'. Make sure it exists"
            " and is executable."
        )
    # Run liftOver
    output = subprocess.run(
        [
            "./liftOver",
            "-positions",
            "liftover_positions.txt",
            chain_file,
            "liftover_output.txt",
            "unMapped",
        ],
        check=True,
        capture_output=True,
    )

    assert output.returncode == 0, (
        "\n\tError in running liftover with file liftover.bed"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )


def merge_in_coords(count_df, liftover_coords_df):
    """
    Merge in the GRCh38 coordinates to the count dataframe

    Parameters
    ----------
    count_df : pd.DataFrame
        dataframe with one row per variant and aggregate counts
    liftover_coords_df : pd.DataFrame
        dataframe with columns for GRCh38 coords

    Returns
    -------
    merged : pd.DataFrame
        dataframe of counts with GRCh38 coordinates added
    """
    merged = pd.concat([count_df, liftover_coords_df], axis=1)

    return merged


def main():
    args = parse_args()
    count_df = read_in_to_df(args.input, sep="\t", header=0)
    create_positions_file_from_count_df(count_df)
    run_liftover(args.chain)
    lifted_over = pd.read_csv(
        "liftover_output.txt", sep="\t", header=None, names=["GRCh38_position"]
    )
    split_lifted_over_coords = split_liftover_bed_columns(
        lifted_over, "GRCh38_position"
    )
    merged = merge_in_coords(count_df, split_lifted_over_coords)
    merged.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
