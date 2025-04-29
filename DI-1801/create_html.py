import argparse
import json
import pandas as pd
from jinja2 import Environment, FileSystemLoader


from utils import read_in_json


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description="Information required to make HTML file"
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Path to input CSV to add in to HTML template",
    )

    parser.add_argument(
        "--columns",
        required=True,
        type=str,
        help="Path to JSON file with column information to include in HTML",
    )

    parser.add_argument(
        "--template",
        required=True,
        type=str,
        help="Name of HTML template file in /templates dir",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Name of output HTML file",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    csv_data = pd.read_csv(
        args.input,
        keep_default_na=False,
    ).to_dict(orient="records")
    columns = read_in_json(args.columns)
    env = Environment(loader=FileSystemLoader("templates"))
    template = env.get_template(args.template)
    output = template.render(csv_data=csv_data, columns=columns)

    with open(args.output, "w", encoding="utf8") as f:
        f.write(output)


if __name__ == "__main__":
    main()
