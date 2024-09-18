import argparse
import dxpy
import json
import numpy as np
import os
import pandas as pd
import requests
import sys
import time

from collections import defaultdict
from dotenv import load_dotenv
from requests.auth import HTTPBasicAuth


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description='Information required to calculate TATs'
    )

    date_format = "%Y-%m-%d"

    parser.add_argument(
        '-s',
        '--start_date',
        required=True,
        type=str,
        help=(
            f"Date to filter projects created after in "
            f"{date_format.replace(r'%', r'%%')} format e.g. 2023-09-18"
        )
    )

    parser.add_argument(
        '-e',
        '--end_date',
        required=True,
        type=str,
        help=(
            f"Date to filter projects created before in "
            f"{date_format.replace(r'%', r'%%')} format e.g. 2024-09-18"
        )
    )

    parser.add_argument(
        '-a',
        '--assay',
        required=True,
        type=str,
        help='Assay type we are auditing, e.g. MYE'
    )

    parser.add_argument(
        '-o',
        '--outfile_name',
        required=True,
        type=str,
        help='Name of the output Excel file'
    )

    return parser.parse_args()


def find_dx_projects(start_date, end_date, assay_type):
    """
    Find 002 projects for an assay between certain dates

    Parameters
    ----------
    start_date : str
        date to filter projects created after e.g. '2024-04-09'
    end_date : str
        date to filter projects created before e.g. '2024-09-04'
    assay_type : str
        assay type we're auditing

    Returns
    -------
    projects : list
        list of dicts, each with info about a DX project
    """
    projects = list(dxpy.find_projects(
        name=f"002*{assay_type}",
        name_mode="glob",
        level="VIEW",
        describe={
            "fields": {
                "id": True,
                "name": True
            }
        },
        created_after=start_date,
        created_before=end_date
    ))

    assert projects, (
        f"No {assay_type} projects found between {start_date} and {end_date}"
    )

    print(f"Found {len(projects)} 002 {assay_type} projects\n")

    return projects


def find_log_file_in_folder(run_name):
    """
    Find run.{run_id}.lane.all.log file in the relevant folder for the run in
    the 001_Staging_Area52 project

    Parameters
    ----------
    run_name: str
        name of the run

    Returns
    -------
    log_file_info : dict
        dict about the log file in the project
    """
    log_files_found = list(
        dxpy.find_data_objects(
            project='project-FpVG0G84X7kzq58g19vF1YJQ',
            folder=f'/{run_name}/runs',
            name="*.lane.all.log",
            name_mode='glob',
            classname='file',
            describe={
                'fields': {
                    'name': True,
                    'created': True
                }
            }
        )
    )

    assert len(log_files_found) == 1, (
        f'No or multiple log files found for run {run_name}'
    )

    log_file_info = log_files_found[0]

    return log_file_info


def get_datetime_log_file_created(log_file_info):
    """
    Convert time file created in DX epoch time to timestamp

    Parameters
    ----------
    log_file_info : dict
        info about the log file for a run

    Returns
    -------
    timestamp_created : str
        timestamp the log file was created
    """
    created_epoch = log_file_info['describe']['created'] / 1000
    timestamp_created = time.strftime(
        '%Y-%m-%d %H:%M:%S', time.localtime(created_epoch)
    )

    return timestamp_created


def query_jira_tickets_in_queue(jira_email, jira_token, queue_id):
    """
    Query helpdesk ticket queue with Jira API. As can't change size of
    response and seems to be limited to 50, loop over each page of response to
    get all tickets until response is empty

    Parameters
    ----------
    jira_email : str
        email for Jira API
    jira_token : str
        token for Jira API
    queue_id :  int
        ID for the relevant servicedesk queue
    Returns
    -------
    response_data :  list
        list of dicts with response from Jira API request
    """
    auth = HTTPBasicAuth(jira_email, jira_token)
    headers = {"Accept": "application/json"}
    base_url = (
        "https://cuhbioinformatics.atlassian.net/rest/servicedeskapi/"
        f"servicedesk/4/queue/{queue_id}/issue"
    )
    response_data = []
    new_data = True
    start = 0
    page_size = 50

    while new_data:
        queue_response = requests.request(
            "GET",
            url=f"{base_url}?start={start}",
            headers=headers,
            auth=auth,
            timeout=600
        )
        # Check request response OK, otherwise exit as would be key error
        if queue_response.ok:
            new_data = json.loads(queue_response.text)['values']
            response_data += new_data
            start += page_size
        else:
            print("Issue with Jira response - check credentials")
            sys.exit(1)

    return response_data


def get_ticket_transition_times(jira_email, jira_token, ticket_id):
    """
    Get times of the Jira ticket transitions to different statuses

    Parameters
    ----------
    ticket_id : int
        the ID of the Jira ticket

    Returns
    -------
    transitions_dict : dict
        dict of the time the ticket transition to a status happened last
        (to account for a ticket being put back to a previous status)
    Example:
    {
        'Data Received': '2024-02-02 09:53:06',
        'Data processed': '2024-02-02 17:09:46',
        'All samples released': '2024-02-05 16:17:08'
    }
    """
    auth = HTTPBasicAuth(jira_email, jira_token)
    headers = {"Accept": "application/json"}

    transitions_dict = defaultdict(list)
    url = (
        "https://cuhbioinformatics.atlassian.net/rest/api/3/issue/"
        f"{ticket_id}/changelog"
    )

    log_response = requests.request(
        "GET",
        url,
        headers=headers,
        auth=auth,
        timeout=600
    )

    change_info = json.loads(log_response.text)['values']

    # Loop over changes, get times the ticket changed to that status
    # append to dict then get the latest time the ticket transitioned
    # to that status
    for change in change_info:
        status_details = [
            x for x in change['items'] if x['field'] == 'status'
        ]

        if status_details:
            status_date, status_time = change['created'].split('T')
            status_time = status_time.split('.', 3)[0]
            status_date_time = f"{status_date} {status_time}"
            new_state = status_details[0]['toString']

            transitions_dict[new_state].append(status_date_time)

    for status_key, status_change_times in transitions_dict.items():
        transitions_dict[status_key] = max(status_change_times)

    return transitions_dict


def add_ticket_resolved_time(
    jira_email, jira_token, jira_api_response, assay_type
):
    """
    Store info about all the tickets we've found for our assay

    Parameters
    ----------
    jira_email : str
        email for Jira API
    jira_token : str
        token for Jira API
    jira_api_response : list
        list of dicts, each with info about a Jira helpdesk ticket
    assay_type : str
        assay type we're auditing

    Returns
    -------
    jira_run_dict : collections.defaultdict
        dict with info about each ticket for our assay
    Example:
    {
        '240916_A01303_0452_BHHWNHDRX5': {
            'ticket_key': 'EBH-3093',
            'ticket_id': '26128',
            'jira_status': 'All samples released',
            'resolved_time': '2024-09-18 10:33:30'
        },
        '240909_A01303_0448_BHJ3VHDRX5': {
            'ticket_key': 'EBH-3068',
            'ticket_id': '26065',
            'jira_status': 'All samples released',
            'resolved_time': '2024-09-12 16:27:24'
        },
        ...
    }
    """
    jira_run_dict = defaultdict(dict)
    for issue in jira_api_response:
        ticket_name  = issue['fields']['summary']
        assay_type_field = issue.get('fields').get('customfield_10070')
        if assay_type_field:
            assay = assay_type_field[0].get('value')
        else:
            assay = 'Unknown'

        # If ticket for our assay, store relevant info about the ticket
        if assay == assay_type:
            jira_run_dict[ticket_name]['ticket_key'] = issue['key']
            jira_run_dict[ticket_name]['ticket_id'] = issue['id']
            jira_run_dict[ticket_name]['jira_status'] = (
                issue['fields']['status']['name']
            )
            change_log = get_ticket_transition_times(
                jira_email, jira_token, issue['id']
            )
            jira_resolved = change_log.get('All samples released')
            if jira_resolved:
                jira_run_dict[ticket_name]['time_resolved'] = jira_resolved

    return jira_run_dict


def create_dict_with_all_info(list_of_projects, ticket_dict):
    """
    Create a final dictionary with all info about a run

    Parameters
    ----------
    list_of_projects : list
        list of dicts, each with info about a DX project
    ticket_dict : collections.defaultdict
        dict containing info about all tickets for the assay we're auditing

    Returns
    -------
    project_dict: collections.defaultdict
        dict with info about upload time and ticket resolved time for each run
    """
    project_dict = defaultdict(dict)

    for project in list_of_projects:
        run_name = project['describe']['name'][4:-4]
        created_epoch = find_log_file_in_folder(run_name)
        created_string = get_datetime_log_file_created(created_epoch)
        project_dict[run_name]['sequence_complete'] = created_string

        if run_name in ticket_dict:
            if ticket_dict[run_name].get('time_resolved'):
                project_dict[run_name]['data_available'] = (
                    ticket_dict[run_name]['time_resolved']
                )
                ticket_key = ticket_dict[run_name]['ticket_key']
                project_dict[run_name]['ticket_link'] = (
                    "https://cuhbioinformatics.atlassian.net/servicedesk/"
                    f"customer/portal/4/{ticket_key}"
                )
            else:
                print(
                    f"Run {run_name} does not have a Jira ticket at "
                    "All samples released"
                )

    return project_dict


def read_to_df(project_dict):
    """
    Read the nested dictionary into a pandas dataframe

    Parameters
    ----------
    project_dict : dict
        defaultdict where each key is a run name and values are info about
        the run

    Returns
    -------
    sorted_run_df: pd.DataFrame
        pandas df with info about each run
    """
    run_df = pd.DataFrame(
        project_dict.values()
    ).assign(run_name=project_dict.keys())

    run_df = run_df[[
        'run_name', 'ticket_link', 'sequence_complete', 'data_available'
    ]]

    # Convert cols to datetime so we can do calculations later
    cols_to_convert = ['sequence_complete', 'data_available']
    run_df[cols_to_convert] = run_df[cols_to_convert].apply(
        pd.to_datetime, format='%Y-%m-%d %H:%M:%S'
    )

    sorted_run_df = run_df.sort_values(
        by='sequence_complete', ignore_index=True
    )

    return sorted_run_df


def add_tat_in_hours(run_df):
    """
    Add the TAT from sequence_complete to data_available for each run in hours

    Parameters
    ----------
    run_df : pd.DataFrame
        pandas df with info about each run

    Returns
    -------
    run_df : pd.DataFrame
        pandas df with info about each run, plus new column 'run_tat'
    """
    run_df['run_tat'] = (
        (run_df['data_available'] - run_df['sequence_complete'])
        / np.timedelta64(1, 'h')
    )

    return run_df

def calculate_mean_and_median_tat(run_df):
    """
    Calculate the mean and median TAT over all runs

    Parameters
    ----------
    run_df : pd.DataFrame
        pandas df with info about tat for each run

    Returns
    -------
    mean_overall_tat : float
        mean TAT of all runs in hours
    median_overall_tat : float
        median TAT of all runs in hours
    """
    mean_overall_tat = run_df['run_tat'].mean()
    median_overall_tat = run_df['run_tat'].median()

    print(
        f"The mean TAT for {len(run_df)} runs is {mean_overall_tat:.2f} hours"
    )
    print(
        f"The median TAT for {len(run_df)} runs is {median_overall_tat:.2f} "
        "hours"
    )

    return mean_overall_tat, median_overall_tat

def write_out_to_excel(tat_df, outfile_name):
    """
    Write out the underlying data to an Excel file

    Parameters
    ----------
    tat_df : pd.DataFrame
        pandas df with info about tat for each run
    outfile_name : str
        name of the Excel file to be written out
    """
    writer = pd.ExcelWriter(outfile_name)
    tat_df.to_excel(writer, sheet_name='Sheet1', index=False)

    # Auto adjust column widths
    for column in tat_df:
        column_length = max(
            tat_df[column].astype(str).map(len).max(), len(column)
        )
        col_idx = tat_df.columns.get_loc(column)
        writer.sheets['Sheet1'].set_column(col_idx, col_idx, column_length)

    writer.save()


def main():
    """Main function"""
    load_dotenv()
    args = parse_args()
    jira_email = os.environ.get('JIRA_EMAIL')
    jira_token = os.environ.get('JIRA_TOKEN')

    # Get all 002 projects for assay within our audit period
    projects_002 = find_dx_projects(
        args.start_date, args.end_date, args.assay
    )

    # Query closed sequencing runs queue to get all Jira tickets
    closed_sequencing_jira_tickets = query_jira_tickets_in_queue(
        jira_email, jira_token, 35
    )
    ticket_dict = add_ticket_resolved_time(
        jira_email, jira_token, closed_sequencing_jira_tickets, args.assay
    )

    # Create final dict with info about each run
    project_dict = create_dict_with_all_info(projects_002, ticket_dict)

    # Convert to df, add TATs and write out to Excel
    run_df = read_to_df(project_dict)
    tat_df = add_tat_in_hours(run_df)
    mean_tat, median_tat = calculate_mean_and_median_tat(tat_df)
    write_out_to_excel(tat_df, args.outfile_name)


if __name__ == "__main__":
    main()
