# This is a temporary placeholder to generate the exons per transcript
# using the zetta url method rather than the API client as it is still
# out of date
import json
import pandas as pd
import urllib.request


def query_cellbasedict(exon_dict):
    txs_dict = {}
    txs_dict['chr'] = 'chr' + exon_dict["chromosome"]
    txs_dict['exon_start'] = exon_dict["genomicCodingStart"]
    txs_dict['exon_end'] = exon_dict["genomicCodingEnd"]
    txs_dict['gene_symbol'] = data['responses'][0]['results'][0]['name']
    # cannot take from exon as it has the _exonnumber attached to the transcript
    txs_dict['transcript_id'] = data['responses'][0]['results'][0]['id']
    txs_dict['exonNumber'] = exon_dict["exonNumber"]
    return txs_dict


transcript_list = ['NM_014915', 'NM_000633',
                'NM_000061', 'NM_003467', 'NM_002755',
                'NM_002661', 'NM_003334', 'NM_014953',
                'NM_017709', 'NM_002460', 'NM_032415', 'NM_006015',
                'NM_001145785', 'NM_002015', 'NM_001664']

# main function
# there's a chance the transcript does not exist in cellbase so let's
# track that
txs_notin_cellbase = []
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
        exons_num = len(data['responses'][0]['results'][0]['exons'])

        # loop through each exon for a transcript and make a long list of dict
        list_txs_dict = []
        for exon in range(exons_num):
            exon_dict = data['responses'][0]['results'][0]['exons'][exon]
            # print(exon_dict)
            if exon_dict["phase"] != -1:
                extracted_exons_dict = query_cellbasedict(exon_dict)
                list_txs_dict.append(extracted_exons_dict)
            else:
                print("exon " + str(exon) + " is a non coding exon")

        # combine all exons into one table for that transcript
        transcript_table = pd.DataFrame(list_txs_dict)

        # append to main df of all transcripts
        df = df.append(transcript_table)

df[['exon_start']] = df[['exon_start']] - 5
df[['exon_end']] = df[['exon_end']] + 5

df.to_csv("coding_unrestricted_GRCh38_myeloid_v2.0_new_capture_regions_only.bed", sep="\t",
        header=False, index=False)

df2 = df[['chr', 'exon_start', 'exon_end', 'transcript_id']]
df2.to_csv("coding_unrestricted_athena_GRCh38_myeloid_v2.0_new_capture_regions_only.bed", sep="\t",
        header=False, index=False)

df.to_csv("exons_nirvana_GRCh38_5bp_flank_v2.0_new_capture_regions_only.tsv", sep="\t",
        header=False, index=False)
