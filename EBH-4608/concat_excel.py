# want to concatonate all of the excel files for haemonc run 002_251121_A01295_0638_BHFWLTDRX7_MYE
import dxpy
import pandas as pd
import openpyxl
import argparse
import subprocess
import re
import io as io

project_id = "project-J4V0y18454QjYbJjzY09XxQp"
file_search_term = ".*xlsx"

# find the necessary files
def find_files_in_project(project: dict, name: str,
                          name_mode: str = "regexp") -> list:
    """
    Find files in a project by name.

    Parameters
    ----------
    project : dict
        The project object.
    name : str
        The name to search for.
    name_mode : str
        The mode for name matching.

    Returns
    -------
    list
        A list of excel file objects matching the name.
    """
    excels = list(
        dxpy.find_data_objects(
            project=project_id,
            name=name,
            name_mode=name_mode,
            classname="file",
            describe={
                "fields": {
                    "archivalState": True,
                    "name": True,
                    "created": True,
                }
            },
        )
    )

    return excels

# convert file information into dataframe
def convert_to_df(
    excel_list: list
) -> pd.DataFrame:
    """
    Convert a list of excel file metadata to a pandas DataFrame, exclude rows
    from DF using exclude_projects and exclude_samples lists.

    Parameters
    ----------
    excel_list : list
        list of excel file metadata dicts.
    exclude_projects : list
        List of project IDs to exclude from the DataFrame.
    exclude_samples : list
        List of sample names or patterns to exclude from the DataFrame.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the excel file metadata.
    """
    flat_list = [
        {
            "file_id": item["id"],
            "project_id": item["project"],
            "project_name": item["project_name"],
            "name": item["describe"]["name"],
            "created": item["describe"]["created"],
            "archive_state": item["describe"]["archivalState"],
        }
        for item in excel_list
    ]
    df = pd.DataFrame(flat_list)

    # Add useful fields for later
    df["sample"] = df["name"].str.split("-").str[0:2].str.join("-")
    df["project_file"] = df["project_id"] + ":" + df["file_id"]
    df["assay"] = df["project_name"].str.split("_").str[-1]

    return df

### Gather excel file IDs ###

excels = []
excels_in_project = find_files_in_project(
            project_id, file_search_term
        )

for excel in excels_in_project:
     excel["project_name"] = "002_251121_A01295_0638_BHFWLTDRX7_MYE"
excels.extend(excels_in_project)

all_excels = convert_to_df(excels)

# save excel file detail 
all_excels.to_csv(
        f"all_run_excels.tsv",
        sep="\t",
        index=False,
    )
print(all_excels)


### combine data drom included and excluded sheets  of each excel ##

# create empty dataframe
merged_excel = pd.DataFrame()

for file_id in all_excels["file_id"]:
    print(file_id)

    # Use dx cat as dxpy.open_dxfile method would lead to
    # hyperlinks being removed for all cells
    cmd = (
       f"dx cat {file_id} "
    )
    output = subprocess.run(cmd, shell=True,
                             capture_output=True, check=False)
    # use openpyxl for hyperlink functionality
    sample_excel = openpyxl.load_workbook(io.BytesIO(output.stdout),
                                           data_only=False)

    # sheets of interest are included and excluded
    for sheet in ["included", "excluded"]:
        sample_sheet = sample_excel[sheet]
        max_row = sample_sheet.max_row
        # remove hyperlink from gnomad (row 20 and 21)
        for row in sample_sheet.iter_cols(min_col=20, 
                                          max_col=21,
                                          min_row = 0,
                                          max_row = max_row):
            #for row in sample_sheet.iter_rows(min_row=0, max_row=max_row):
            for cell in row:
                # cell.hyperlink showing none and cell.style showing normal
                # cell.value showing a hyperlink string
                if "=HYPERLINK"in str(cell.value):
                    print(cell.value)
                    # extract value from hyperlink strin
                    new_value = str(cell.value).split(",")[1].split(")")[0]
                    cell.value = float(new_value)
                    print(cell.value)
        # convert to dataframe 
        sample_data = sample_sheet.values
        colnames = [c.value for c in next(sample_sheet.iter_rows(min_row=1,
                                                                 max_row=1))]
        sample_tab = pd.DataFrame(sample_data, columns=colnames)
        # remove first row because duplicated index
        sample_tab = sample_tab.drop([0])
        # remove where sampleName is null as pyexcel also loads empty rows
        sample_tab = sample_tab[~sample_tab["sampleName"].isnull()]
        sample_tab["Included_or_Excluded"] = sheet
        merged_excel =pd.concat([merged_excel,sample_tab])

# move included or excluded column to front so easily visible
included_excluded_col = merged_excel.pop("Included_or_Excluded")
merged_excel.insert(1, "Included_or_Excluded", included_excluded_col)

# save to excel
print(merged_excel)
merged_excel.to_excel("002_251121_A01295_0638_BHFWLTDRX7_MYE_merged_workbooks.xlsx",
                       index=False)

