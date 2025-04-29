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
        description="Information required to normalise variants"
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
        help="Path to JSON file with info fields for VCF",
    )

    parser.add_argument(
        "--output_csv",
        required=True,
        type=str,
        help="Name of output CSV file with aggregated counts",
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
    Read in fASTA to pysam.FastaFile object

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


def get_ref_contig_info(fasta):
    """
    Get a list of contigs from a FASTA file

    Parameters
    ----------
    fasta : pysam.FastaFile
        FASTA file as pysam object

    Returns
    -------
    contigs : list
        list of strings, each representing a FASTA contig
    """
    contigs = []
    for contig in fasta.references:
        length = fasta.get_reference_length(contig)
        contigs.append(f"##contig=<ID={contig},length={length}>")
    return contigs


def normalise_variant_row(row, fasta):
    """
    Normalise a row of data

    Parameters
    ----------
    row : pd.Series
        Row of data to normalise
    fasta : pysam.FastaFile
        FASTA file as pysam object

    Returns
    -------
    norm_chrom : str
        Normalised chrom
    norm_pos : int
        Normalised position
    norm_ref : str
        Normalised reference allele
    norm_alt : str
        Normalised alternate allele
    """
    # Set explicit datatypes to avoid any issues when querying
    chrom, pos, ref, alt = (
        str(row["Chromosome"]),
        int(row["Start_Position"]),
        str(row["Reference_Allele"]),
        str(row["Tumor_Seq_Allele2"]),
    )

    # Set norm values to original by default
    norm_chrom, norm_pos, norm_ref, norm_alt = chrom, pos, ref, alt

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
            norm_ref = prefix_bp
            norm_pos = pos
            norm_alt = prefix_bp + alt
        # It's a deletion - we need the remove one from position
        else:
            norm_ref = prefix_bp + ref
            norm_alt = prefix_bp
            norm_pos = pos - 1

    # For non-indels, verify ref allele
    else:
        ref_from_fasta = ref_seq[1 : len(ref) + 1]
        if ref != ref_from_fasta:
            # Reference mismatch - keep original values
            print(
                f"Reference mismatch for {chrom}-{pos}-{ref}-{alt} -"
                f" {ref_from_fasta} found in FASTA"
            )

    return norm_chrom, norm_pos, norm_ref, norm_alt


def convert_to_vcf_representation(genie_count_data, fasta):
    """
    Make new columns in our df with the normalised chrom, pos, ref and alt in
    GRCh37

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
        chrom_norm, pos_norm, ref_norm and alt_norm
    """
    genie_count_data[["chrom_norm", "start_norm", "ref_norm", "alt_norm"]] = (
        genie_count_data.apply(
            lambda row: pd.Series(normalise_variant_row(row, fasta)), axis=1
        )
    )

    return genie_count_data


def write_variants_to_vcf(
    genie_counts_normalised, output_vcf, fasta, info_fields
):
    """
    _summary_

    Parameters
    ----------
    genie_counts_normalised : _type_
        _description_
    output_vcf : _type_
        _description_
    fasta : _type_
        _description_
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
        "Genie_original_coords",
        "1",
        "String",
        "Original variant info from Genie",
    )

    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)
    for _, row in genie_counts_normalised.iterrows():
        # Create a new record
        record = vcf_out.new_record(
            contig=str(row["chrom_norm"]),
            start=int(row["start_norm"]) - 1,
            alleles=(row["ref_norm"], row["alt_norm"]),
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
                    field_value = field_value.replace(",", "&")
                record.info[field["id"]] = field_value
        # Add in original Genie GRCh37 chrom-pos-ref-alt
        orig_coord_str = f"{row['Chromosome']}_{row['Start_Position']}_{row['Reference_Allele']}_{row['Tumor_Seq_Allele2']}"
        record.info["Genie_original_coords"] = orig_coord_str
        # Write the record
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
    genie_counts_normalised = convert_to_vcf_representation(
        genie_counts, fasta
    )
    genie_counts_normalised.to_csv(
        args.output_csv,
        index=False,
    )
    info_fields = read_in_json(args.info_fields)
    write_variants_to_vcf(
        genie_counts_normalised, args.output_vcf, fasta, info_fields
    )


if __name__ == "__main__":
    main()
