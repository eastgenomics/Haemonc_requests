# This is a temporary placeholder to generate the exons per transcript
# using the zetta url method rather than the API client as it is still
# out of date
import argparse
import json
import pandas as pd
import urllib.request

def parse_args():
    """Parse through input arguments
    Returns:
        args (Namespace): Namespace that you can extract relevant input arguments
    """
    # Read in arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-t', '--transcript_file',
        help='Transcripts listed in a tsv',
        required=True
        )

    parser.add_argument(
        '-v', '--version_number',
        type=int,
        help='What version of the HaemOnc transcript is this? (integer)',
        required=True
        )

    args = parser.parse_args()

    return args


def query_cellbasedict(exon_dict, data):
    """Queries info from exon dict

    Args:
        exon_dict (dictionary): dictionary containing all info
        about the exon. Usually take the coding regions of exons.
        data (dictionary): dictionary containing all informatio about
        the transcript.

    Returns:
        txs_dict (dictionary): a dictionary of selected exon info
    """
    txs_dict = {}
    txs_dict['chr'] = 'chr' + exon_dict["chromosome"]
    txs_dict['exon_start'] = exon_dict["genomicCodingStart"]
    # We need to include the UTR region (negative strand so its the end)
    # of the ANKRD26 gene. We will use the exon[end] rather than
    # exon[genomicCodingEnd]
    if data['responses'][0]['results'][0]['name'] == "ANKRD26" and exon_dict["exonNumber"] == 1:
        txs_dict['exon_end'] = exon_dict["end"]
    else:
        txs_dict['exon_end'] = exon_dict["genomicCodingEnd"]
    txs_dict['gene_symbol'] = data['responses'][0]['results'][0]['name']
    # cannot take from exon as it has the _exonnumber attached to the transcript
    txs_dict['transcript_id'] = data['responses'][0]['results'][0]['id']
    txs_dict['exonNumber'] = exon_dict["exonNumber"]
    return txs_dict

def main():
    """Generates a file of chr, start, end, gene symbol,
       transcript and exon number from a list of transcripts.

       Outputs three files, where coding_unrestricted_GRCh38_myeloid and
       exons files are the same but these files are appended to other
       files so it's seperate outputs makes it easier for appending files.
    """
    args = parse_args()
    with open(args.transcript_file, "r") as txs_file:
        transcript_list = txs_file.read().splitlines()
        transcript_list = list(filter(None, transcript_list))

    df = pd.DataFrame()

    for transcript in transcript_list:
        # get transcript from cellbase via webserver
        print(transcript)
        url_get = "https://ws.zettagenomics.com/cellbase/webservices/rest/v5/hsapiens/feature/transcript/"
        url_end = "/info?source=refseq"
        url_address = url_get + transcript + url_end
        with urllib.request.urlopen(url_address) as url:
            data = json.loads(url.read().decode())

        if not data['responses'][0]['results']:
            print("RefSeq transcript does not exist in cellbase")
        else:
            # number of exons:
            all_exons = data['responses'][0]['results'][0]['exons']
            # loop through each exon for a transcript and make a long list of dict
            list_txs_dict = []
            for exon in all_exons:
                if exon["phase"] != -1:
                    # is the exon phase is -1, this means the whole exon
                    # is not translated so we can skip these. also the
                    # start and end will be 0, so its easier to skip.
                    extracted_exons_dict = query_cellbasedict(exon, data)
                    list_txs_dict.append(extracted_exons_dict)

            # combine all exons into one table for that transcript
            transcript_table = pd.DataFrame(list_txs_dict)

            # append to main df of all transcripts
            df = df.append(transcript_table)

    df[['exon_start']] = df[['exon_start']] - 5
    df[['exon_end']] = df[['exon_end']] + 5

    df.to_csv("coding_unrestricted_GRCh38_myeloid_5bp_flank_v" + str(args.version_number) +".0.0.bed", sep="\t",
            header=False, index=False)

    df2 = df[['chr', 'exon_start', 'exon_end', 'transcript_id']]
    df2.to_csv("coding_unrestricted_athena_GRCh38_myeloid_5bp_flank_v" + str(args.version_number) +".0.0.bed", sep="\t",
            header=False, index=False)

    df.to_csv("exons_cellbase_GRCh38_5bp_flank_v" + str(args.version_number) +".0.0.bed", sep="\t",
            header=False, index=False)

if __name__ == "__main__":

    main()