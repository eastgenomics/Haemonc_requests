import argparse
from Bio import Entrez
import pandas as pd

from collections import defaultdict
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
        description="Information required to subset MAF file"
    )

    parser.add_argument(
        "--input_maf",
        required=True,
        type=str,
        help="Path to MAF file unmodified from Genie",
    )

    parser.add_argument(
        "--bed_file",
        required=True,
        type=str,
        help=(
            "Path to haemonc BED file which has gene symbols in the 5th column"
        ),
    )

    parser.add_argument(
        "--entrez_email",
        required=True,
        type=str,
        help="Email address to identify yourself for Entrez queries",
    )

    parser.add_argument(
        "--extra_symbols",
        required=False,
        type=str,
        help=(
            "Extra gene symbols to include (comma separated) e.g. 'TP53,EGFR'"
        ),
    )

    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Name of output MAF file subsetted to haemonc genes",
    )

    return parser.parse_args()


def get_entrez_gene_id(gene_symbol):
    """
    Get a list of Entrez gene ID(s) for a gene symbol

    Parameters
    ----------
    gene_symbol : str
        gene symbol to query

    Returns
    -------
    list
        list of Entrez IDs (list of strings) or None if none found
    """
    handle = Entrez.esearch(
        db="gene", term=f"{gene_symbol}[Sym] AND Homo sapiens[Organism]"
    )
    record = Entrez.read(handle)
    handle.close()
    if record["IdList"]:
        if len(record["IdList"]) > 1:
            print(
                f"Gene {gene_symbol} has more than one Entrez ID, searching"
                " with all"
            )
            print(record["IdList"])
        return record["IdList"]
    return None


def get_haemonc_gene_ids(gene_symbols):
    """
    Query a list of gene symbols to get all their Entrez gene IDs

    Parameters
    ----------
    gene_symbols : list
        list of strings, each a gene symbol to query

    Returns
    -------
    entrez_ids : collections.defaultdict(list)
        dict with each gene and its Entrez gene IDs as value (list of ints)
    """
    entrez_ids = defaultdict(list)
    for gene in gene_symbols:
        entrez_id = get_entrez_gene_id(gene)
        if not entrez_id:
            print(f"No Entrez gene ID found for gene {gene}")
        else:
            entrez_ids[gene] = list(map(int, entrez_id))

    return entrez_ids


def check_gene_matching(haemonc_gene_symbols, genie_data):
    """
    Check what has been kept in the Genie data, either by gene symbol or
    Entrez ID

    Parameters
    ----------
    haemonc_gene_symbols : collections.defaultdict(list)
        dict where key is gene symbol and value is list of entrez IDs
    genie_data : pd.DataFrame
        dataframe of all genie data (unmodified)
    """
    # Make reverse mapping of Entrez IDs to gene symbols
    entrez_to_symbols = defaultdict(set)
    for symbol, ids in haemonc_gene_symbols.items():
        for eid in ids:
            entrez_to_symbols[eid].add(symbol)

    all_symbols = set(haemonc_gene_symbols.keys())
    all_entrez_ids = set(entrez_to_symbols.keys())

    # Identify which symbols/IDs are matched in the Genie dataset
    observed_symbols = set(genie_data["Hugo_Symbol"])
    observed_ids = set(genie_data["Entrez_Gene_Id"].dropna())

    # Identify gene symbols common to both datasets
    matched_by_symbol = all_symbols & observed_symbols
    # Find Haemonc gene symbols that match based on ID
    matched_by_id = {
        symbol
        for symbol, ids in haemonc_gene_symbols.items()
        if any(eid in observed_ids for eid in ids)
    }
    # Find symbols not matched based on symbol or Entrez ID
    matched_neither = all_symbols - (matched_by_symbol | matched_by_id)

    # Print genes not matched at all
    if matched_neither:
        print("\nDid not match by either Hugo_Symbol or Entrez_Gene_Id:")
        for symbol in sorted(matched_neither):
            print(f"{symbol}: {haemonc_gene_symbols[symbol]}")

    # Mismatches: ID matches, but symbol is unexpected
    mask_id_match = genie_data["Entrez_Gene_Id"].isin(all_entrez_ids)
    mask_symbol_unexpected = ~genie_data["Hugo_Symbol"].isin(all_symbols)
    mismatched = genie_data[mask_id_match & mask_symbol_unexpected].copy()

    mismatched["Expected_Symbols"] = mismatched["Entrez_Gene_Id"].map(
        lambda eid: ", ".join(sorted(entrez_to_symbols.get(eid, [])))
    )
    unique_mismatches = mismatched.drop_duplicates(
        subset=["Hugo_Symbol", "Entrez_Gene_Id"]
    )

    if not unique_mismatches.empty:
        print(
            "\nMismatched rows (Entrez_Gene_Id belongs to a different symbol):"
        )
        for row in unique_mismatches.itertuples(index=False):
            print(
                f"Hugo_Symbol = '{row.Hugo_Symbol}', Entrez_Gene_Id ="
                f" {row.Entrez_Gene_Id} => Expected symbol(s):"
                f" '{row.Expected_Symbols}'"
            )


def subset_data(haemonc_entrez_ids, genie_data):
    """
    Subset the Genie data by the haemonc gene symbols and Entrez IDs

    Parameters
    ----------
    haemonc_entrez_ids : dict
        dict of each gene symbol and its Entrez IDs
    genie_data : pd.DataFrame
        Genie data to subset

    Returns
    -------
    genie_haemonc_genes : pd.DataFrame
        subset of Genie data for just haemonc genes
    """
    # Create set of all Entrez IDs for our symbols
    entrez_ids = {eid for ids in haemonc_entrez_ids.values() for eid in ids}
    my_gene_symbols = set(haemonc_entrez_ids.keys())
    # Find rows with matching Entrez IDs
    matched_by_id = genie_data[genie_data["Entrez_Gene_Id"].isin(entrez_ids)]

    # Collect Hugo_Symbols from those rows (because some rows with the same
    # gene symbol can weirdly have NAs in the Entrez ID column)
    symbols_from_id_match = set(matched_by_id["Hugo_Symbol"].dropna())

    # Combine symbols above with our original list
    all_symbols_to_keep = my_gene_symbols.union(symbols_from_id_match)

    genie_haemonc_genes = genie_data[
        genie_data["Hugo_Symbol"].isin(all_symbols_to_keep)
    ].reset_index(drop=True)

    return genie_haemonc_genes


def main():
    args = parse_args()
    haemonc_bed = read_in_to_df(args.bed_file, header=None)
    # Get list of unique gene symbols from the BED file
    haemonc_gene_symbols = list(set(haemonc_bed[4].tolist()))
    print(f"{len(haemonc_gene_symbols)} haemonc genes (by symbol) in BED file")
    if args.extra_symbols:
        extra_symbols = [
            item.strip() for item in args.extra_symbols.split(",")
        ]
        print(
            f"{len(extra_symbols)} extra gene symbol(s) to include:"
            f" {extra_symbols}\n"
        )
        haemonc_gene_symbols = haemonc_gene_symbols + extra_symbols
    Entrez.email = args.entrez_email
    haemonc_entrez_ids = get_haemonc_gene_ids(haemonc_gene_symbols)

    genie_data = read_in_to_df(
        args.input_maf,
        header=0,
    )
    # Convert Entrez Gene ID, if exists, to int
    genie_data["Entrez_Gene_Id"] = pd.to_numeric(
        genie_data["Entrez_Gene_Id"], errors="coerce"
    ).astype("Int64")

    check_gene_matching(haemonc_entrez_ids, genie_data)

    genie_haemonc_genes = subset_data(haemonc_entrez_ids, genie_data)
    genie_genes_by_symbol_and_id = list(
        genie_haemonc_genes["Hugo_Symbol"].unique()
    )
    print(
        f"\n{len(genie_genes_by_symbol_and_id)} unique gene symbols in"
        f" subsetted Genie data and {len(genie_haemonc_genes)} variant"
        " lines remaining in haemonc genes"
    )

    genie_haemonc_genes.to_csv(
        args.output,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
