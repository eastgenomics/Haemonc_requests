import argparse
import pandas as pd


def parse_args():
    """Reads on argument passed at the cmd line
    Returns:
        args: args held in an object
    """

    parser = argparse.ArgumentParser(
        description='Anonymising tsv file from uranus handler app')
    parser.add_argument(
        '-f', '--file',
        required=True,
        help='tsv file'
        )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='output file name'
        )
    parser.add_argument(
        '-p', '--output_path',
        required=True,
        help='output path name'
        )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    file = args.file
    output = args.output
    output_path = args.output_path

    df = pd.read_csv(file, sep='\t')

    df2 = df.drop(df.columns[0], axis=1)

    output_filename = output_path + "/" + output + ".tsv"

    df2.to_csv(output_filename, sep="\t", header=True, index=False)


if __name__ == "__main__":

    main()
