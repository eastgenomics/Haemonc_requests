import pandas as pd
import json
import sys


def read_in_to_df(
    filename: str,
    sep: str = "\t",
    skiprows: int = 0,
    usecols: list = None,
    header: int = None,
    dtype: dict = None,
) -> pd.DataFrame:
    """
    Read in file to pandas dataframe

    Parameters
    ----------
    filename : str
        Name of file to read in with pandas
    sep : str, optional
        the separator used in the file, by default "\t"
    skiprows : int, optional
        number of rows to skip when reading in, by default 0
    usecols : list, optional
        list of column names to only read in, by default None
    header : int, optional
        row number to use as column names, by default None
    dtype: dict, optional
        dict to specify datatypes of specific columns, by default None

    Returns
    -------
    pd.DataFrame
        dataframe of input data
    """
    try:
        df = pd.read_csv(
            filename,
            sep=sep,
            skiprows=skiprows,
            usecols=usecols,
            header=header,
            dtype=dtype,
        )
        return df
    except Exception as e:
        print(f"Error reading in file: {e}")
        sys.exit(1)


def read_txt_file_to_list(filename: str) -> list:
    """
    Read TXT file to list

    Parameters
    ----------
    filename : str
        name of file to read in to list

    Returns
    -------
    contents_list : list
        list where each line in the file is an entry
    """
    with open(filename, mode="r", encoding="utf8") as f:
        file_contents = f.read().splitlines()

    # Remove any empty lines, leading and trailing whitespace and duplicates
    contents_list = list(
        set(
            [item for item in (line.strip() for line in file_contents) if item]
        )
    )

    return contents_list


def read_in_json(json_file: str) -> list:
    """
    Read in JSON file to list of dictionaries

    Parameters
    ----------
    json_file : str
        Name of JSON file to read in

    Returns
    -------
    data : list
        List of dictionaries from JSON file
    """
    with open(json_file, "r", encoding="utf8") as f:
        data = json.load(f)
    return data
