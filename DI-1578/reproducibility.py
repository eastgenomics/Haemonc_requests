#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:39:39 2024

@author: arun
"""

import dxpy


def create_folder(project_id, folder_path):
    """
    Create a folder for me in project of interest

    Parameters
    ----------
    project_id : str
        DNAnexus project_id
    folder_path : str
        DNAnexus folder path, should start with /

    Returns
    -------
    None.

    """
    dxpy.api.project_new_folder(project_id, {"folder": folder_path,
                                             "parents":True})
    

def run_sentieon_bwa(fastqc_dict, project_id, folder_path, alias):
    """
    Run the sentieon_bwa command for me

    Parameters
    ----------
    fastqc_dict : dict
        Dictionary containing fastqc reads, for example:
            oncospan_fastqc = {"Name":"oncospan-cell-line-1st",
                               "fastq_L001_R1":"file-id",
                               "fastq_L002_R1":"file-id",
                               "fastq_L001_R2":"file-id",
                               "fastq_L002_R2":"file-id"}
    project_id : str
        DNAnexus project ID where you want to store job output
    folder_path : str
        DNAnexus folder path where you want to store job output, start with /
    alias : str
        version number of the sentieon_bwa app

    Returns
    -------
    str
        App will run and the function will return the associated job_id

    """
    
    # Set the job input
    job_input = {
            "reads_fastqgzs": [dxpy.dxlink(fastqc_dict['fastq_L001_R1']),
                               dxpy.dxlink(fastqc_dict['fastq_L002_R1'])],
            "reads2_fastqgzs": [dxpy.dxlink(fastqc_dict['fastq_L001_R2']),
                                dxpy.dxlink(fastqc_dict['fastq_L002_R2'])],
            "genomebwaindex_targz": dxpy.dxlink("file-Gb76f204XGybZ3J6F731xkBp"),
            "genome_fastagz": dxpy.dxlink("file-Gb757784XGyY3FPvkPQ74K9z")
            }
    
    # Ensure folder path is present 
    create_folder(project_id, folder_path)
    
    # Run the job
    job_run = dxpy.DXApp(name="sentieon-bwa", alias=alias).run(
        job_input,
        project=project_id,
        folder=folder_path,
        name=f"sentieon-bwa_v{alias}_{fastqc_dict['Name']}"
        )
    
    # Return the job_id
    return job_run.describe(fields={"id":True})["id"]


def run_sentieon_tnbam_from_job(job_id, folder_path, alias):
    """
    Run sentieon_tnbam from sentieon_bwa_men job_id outputs
    
    Parameters
    ----------
    job_id : str
        Job_id string of sentieon_bwa_men which we want to take its output as 
        input.
    folder_path : str
        DNAnexus folder path where you want to store job output, project set 
        by job origin, string should start with /
    alias : str
        version number of the sentieon_bwa app

    Returns
    -------
    str
        App will run and the function will return the associated job_id

    """
    
    # wait for the job is done in the first place
    job = dxpy.bindings.DXJob(job_id)
    job.wait_on_done()
    
    # set the job input
    bam_file_id = job.describe()['output']['mappings_bam']['$dnanexus_link']
    job_input = {
        'gatk_resource_bundle': 
            dxpy.dxlink('file-F3zvKp84fXX8qJx43zZXP395',
                      project_id='project-F3zqGV04fXX5j7566869fjFq'),
        'genome_fastagz': dxpy.dxlink('file-Gb757784XGyY3FPvkPQ74K9z'),
        'run_bqsr': True,
        'somatic_algo': 'TNhaplotyper2',
        'mappings_bam': dxpy.dxlink(bam_file_id)
        }
    
    # Ensure folder path is present
    project_id = job.describe()['project']
    create_folder(project_id, folder_path)
    
    # Run the job
    job_name = f"sentieon-tnbam_v{alias}_from_{job.describe()['name']}"
    job_run = dxpy.DXApp(name="sentieon-tnbam", alias=alias).run(
        job_input,
        project=project_id,
        folder=folder_path,
        name=job_name)
    
    # return the job_id
    return job_run.describe(fields={"id":True})["id"]
    
    

def main():
    # Set the main oncospan fastqc for the test
    oncospan_fastqc = {"Name":"oncospan-cell-line-1st",
                       "fastq_L001_R1":"file-G2X22Jj4qyfqzpybKfYJ15VQ",
                       "fastq_L002_R1":"file-G2X25Q04qyfz1qf57gjxqx5f",
                       "fastq_L001_R2":"file-G2X22KQ4qyfgpvyy9GQBkgF9",
                       "fastq_L002_R2":"file-G2X25Qj4qyfZ77jf5kkYV88v"}
    
    # set standard variables
    standard_path = "/reproducibility"
    standard_project = "project-Gjz0P404fV6bq4J30qBfPBFf"
    
    # run sentieon bwa 10 times and store jobs in a list
    bwa_mem_job_ids = []
    for i in range(0,10):
        job_id = run_sentieon_bwa(oncospan_fastqc,
                                  project_id= standard_project,
                                  folder_path=standard_path,
                                  alias="4.2.2")
        bwa_mem_job_ids.append(job_id)
    
    print("All bwa-mem jobs have started running")
    # run sentieon tnbam for each bwa-mem job-id, sthore new jobs in a list
    tnbam_job_ids = []
    for bwa_mem_job in bwa_mem_job_ids:
        tnbam_job_id = run_sentieon_tnbam_from_job(bwa_mem_job,
                                                   folder_path=standard_path,
                                                   alias="5.0.1")
        tnbam_job_ids.append(tnbam_job_id)

    # Print me out the jobs that I want
    print("Here are the list of jobs for bwa-mem:")
    print(bwa_mem_job_ids)
    print("\n")
    print("Here are the list of jobs for tnbam:")
    print(tnbam_job_ids)


if __name__=="__main__":
    main()
    


