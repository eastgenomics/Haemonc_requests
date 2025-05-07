import argparse
import os
import pandas as pd
import pysam
import re

from utils import read_in_to_df, read_in_json


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description=(
            "Information required to convert variants to VCF description"
            " and write out as CSV and VCF files"
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Path to file with variant counts (one row per variant)",
    )

    parser.add_argument(
        "--fasta",
        required=True,
        type=str,
        help="Path to FASTA file used to get ref alleles",
    )

    parser.add_argument(
        "--info_fields",
        required=True,
        type=str,
        help=(
            "Path to JSON file with columns to be added as INFO fields in VCF"
        ),
    )

    parser.add_argument(
        "--output_csv",
        required=True,
        type=str,
        help=(
            "Name of output CSV file with aggregated counts and VCF-like"
            " description"
        ),
    )

    parser.add_argument(
        "--output_vcf",
        required=True,
        type=str,
        help="Name of output VCF file",
    )

    return parser.parse_args()


def read_in_fasta(filename):
    """
    Read in FASTA to pysam.FastaFile object

    Parameters
    ----------
    filename : str
        path to FASTA file

    Returns
    -------
    fasta : pysam.FastaFile
        FASTA file as pysam object
    Raises
    ------
    FileNotFoundError
        If FASTA file does not exist
    ValueError
        If FASTA file is not in the correct format
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"FASTA file not found: {filename}")
    try:
        fasta = pysam.FastaFile(filename)
        return fasta
    except ValueError as err:
        raise ValueError(f"Invalid FASTA format: {err}") from err


def convert_maf_like_variant_to_vcf_description(row, fasta):
    """
    Convert a row of data from MAF-like format to VCF-like format

    Parameters
    ----------
    row : pd.Series
        Row of data to convert
    fasta : pysam.FastaFile
        FASTA file as pysam object

    Returns
    -------
    vcf_chrom : str
        Chrom in VCF-like format
    vcf_pos : int
        Position in VCF-like format
    vcf_ref : str
        Reference allele in VCF-like format
    vcf_alt : str
        Alt allele in VCF-like format
    """
    # Set explicit datatypes to avoid any issues when querying
    chrom, pos, ref, alt = (
        str(row["Chromosome"]),
        int(row["Start_Position"]),
        str(row["Reference_Allele"]),
        str(row["Tumor_Seq_Allele2"]),
    )

    # Set VCF-like values to original by default
    vcf_chrom, vcf_pos, vcf_ref, vcf_alt = chrom, pos, ref, alt

    # Replace ref and alt '-' or other weird chars with ""
    ref = "" if re.fullmatch(r"[-?0]+", ref) else ref
    alt = "" if re.fullmatch(r"[-?0]+", alt) else alt

    # MAF coords are 1-based, pysam uses 0-based half-open intervals and
    # for deletions we want the base before the start, so needs pos - 2
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.FastaFile.fetch
    ref_seq = ""
    ref_seq = fasta.fetch(chrom, pos - 2, pos + len(ref)).upper()
    if not ref_seq:
        print(
            "No reference sequence found in FASTA for"
            f" {chrom}-{pos}-{ref}-{alt}"
        )

    # It's an indel
    if (len(ref) == 0) or (len(alt) == 0):
        prefix_bp = ref_seq[0]
        # For simple insertions, pos is already position of preceding bp
        if ref == "" and len(ref_seq) > 0:
            prefix_bp = ref_seq[1]
            vcf_ref = prefix_bp
            vcf_pos = pos
            vcf_alt = prefix_bp + alt
        # It's a deletion - we need to remove 1 from position
        else:
            vcf_ref = prefix_bp + ref
            vcf_alt = prefix_bp
            vcf_pos = pos - 1

    # For non-indels, verify ref allele
    else:
        ref_from_fasta = ref_seq[1 : len(ref) + 1]
        if ref != ref_from_fasta:
            print(
                f"Reference mismatch for {chrom}-{pos}-{ref}-{alt} -"
                f" {ref_from_fasta} found in FASTA"
            )

    return vcf_chrom, vcf_pos, vcf_ref, vcf_alt


def convert_to_vcf_representation(genie_count_data, fasta):
    """
    Make new columns with the VCF-like description for chrom, pos,
    ref and alt in GRCh37

    Parameters
    ----------
    genie_count_data : pd.DataFrame
        Genie data with one row per variant and count info
    fasta : pysam.FastaFile
        FASTA file as pysam object

    Returns
    -------
    genie_count_data : pd.DataFrame
        Genie data with one row per variant and count info, with new columns
        chrom_vcf, pos_vcf, ref_vcf and alt_vcf
    """
    genie_count_data[["chrom_vcf", "pos_vcf", "ref_vcf", "alt_vcf"]] = (
        genie_count_data.apply(
            lambda row: pd.Series(
                convert_maf_like_variant_to_vcf_description(row, fasta)
            ),
            axis=1,
        )
    )

    return genie_count_data


def write_variants_to_vcf(
    genie_counts_vcf_description, output_vcf, fasta, info_fields
):
    """
    Write out the dataframe (with VCF-like description) to a VCF file
    with the INFO fields specified in the JSON file

    Parameters
    ----------
    genie_counts_vcf_description : pd.DataFrame
        Dataframe with one row per variant and count info and columns with
        VCF-like description
    output_vcf : str
        Name of output VCF file
    fasta : pysam.FastaFile
        FASTA file as pysam object
    info_fields : list
        List of dictionaries with INFO field names and descriptions
        to be added to the VCF file
    """
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    for contig in fasta.references:
        header.add_line(
            f"##contig=<ID={contig},length={fasta.get_reference_length(contig)}>"
        )
    header.formats.add("GT", "1", "String", "Genotype")
    for field in info_fields:
        header.info.add(
            field["id"], field["number"], field["type"], field["description"]
        )
    header.info.add(
        "Genie_original_description",
        "1",
        "String",
        "Original variant info from Genie",
    )

    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)
    # For each original variant, write new variant record with all INFO
    # fields specified in the JSON file
    # Take 1 away from start due to differences in Pysam representation
    for _, row in genie_counts_vcf_description.iterrows():
        record = vcf_out.new_record(
            contig=str(row["chrom_vcf"]),
            start=int(row["pos_vcf"]) - 1,
            alleles=(row["ref_vcf"], row["alt_vcf"]),
            id=".",
            qual=None,
            filter=None,
        )
        for field in info_fields:
            field_name = field["column_name"]
            if field_name in row and pd.notna(row[field_name]):
                field_type = field["type"]
                if field_type == "String":
                    field_value = str(row[field_name])
                elif field_type == "Integer":
                    field_value = int(row[field_name])

                if field_name == "Consequence":
                    # Replace any commas with "&" to avoid VCF parsing issues
                    # of commas
                    field_value = field_value.replace(",", "&")
                record.info[field["id"]] = field_value

        # Add in original Genie GRCh37 chrom-pos-ref-alt
        orig_coord_str = f"{row['Chromosome']}_{row['Start_Position']}_{row['Reference_Allele']}_{row['Tumor_Seq_Allele2']}"
        record.info["Genie_original_description"] = orig_coord_str

        vcf_out.write(record)

    vcf_out.close()
    fasta.close()


def main():
    args = parse_args()
    genie_counts = read_in_to_df(
        args.input,
        sep=",",
        header=0,
        dtype={
            "Entrez_Gene_Id": "Int64",
            "Start_Position": "Int64",
        },
    )
    fasta = read_in_fasta(args.fasta)
    genie_counts_vcf_desc = convert_to_vcf_representation(genie_counts, fasta)
    genie_counts_vcf_desc.to_csv(
        args.output_csv,
        index=False,
    )
    info_fields = read_in_json(args.info_fields)
    write_variants_to_vcf(
        genie_counts_vcf_desc, args.output_vcf, fasta, info_fields
    )


if __name__ == "__main__":
    main()
