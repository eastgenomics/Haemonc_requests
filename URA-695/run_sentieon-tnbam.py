import dxpy as dx
import json

def get_sample_bam_file(sample):
    """
    Find the .bam files from each run which is an
    output of MultiQC

    Parameters
    ----------
    sample : str
        name of the sample that we are interested in acquiring the bam file

    Returns
    -------
    sample_dict : dict
        dict of bam files for versions v3.2.0 v4.2.2
    """
    # Get objects from DNAnexus
    prod_file_response = dx.find_data_objects(
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        name="^(TMv2OCT20-)?"+ sample +"_S[0-9]{2}_L001_markdup.bam$",
        classname='file',
        folder="/sentieon-bwa_v3.2.0",
        name_mode='regexp')
    
    dev_file_response = dx.find_data_objects(
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        name="^(TMv2OCT20-)?"+ sample +"_S[0-9]{2}_L001_markdup.bam$",
        classname='file',
        folder="/sentieon-bwa_v4.2.2",
        name_mode='regexp')
    
    sample_dict = {}
    sample_dict['Name'] = sample
    
    print(f"Processing {sample}")
    for x in prod_file_response:
        sample_dict['sentieon-bwa_v3.2.0'] = x['id']
        print(f"v3.2.0: {x['id']}")
    
    for x in dev_file_response:
        sample_dict['sentieon-bwa_v4.2.2'] = x['id']
        print(f"v4.2.2: {x['id']}")

    return sample_dict


sample_list = ["oncospan-cell-line-1st", "oncospan-cell-line-2nd", 
               "HD734-1st", "HD734-2nd", "HD829-1st", "HD829-2nd", 
               "129404211-24107K0073-24NGSHO18-8128-M-96527893",
               "129459591-24109K0055-24NGSHO18-8128-F-96527893",
               "129450697-24109K0025-24NGSHO18-8128-M-96527893",
               "129325254-24102K0083-24NGSHO17-8128-F-96527893"]

bamfile_list =[]
for sample in sample_list:
    sample_dict = get_sample_bam_file(sample)
    bamfile_list.append(sample_dict)

    
# Make a for loop to make the jobs run in the command line
set_job_ids = []
for sample in bamfile_list:
    job_ids = {}
    job_ids['Name'] = sample['Name']
    
    job_input = {
        'gatk_resource_bundle': 
            dx.dxlink('file-F3zvKp84fXX8qJx43zZXP395',
                      project_id='project-F3zqGV04fXX5j7566869fjFq'),
        'genome_fastagz': dx.dxlink('file-Gb757784XGyY3FPvkPQ74K9z'),
        'run_bqsr': True,
        'somatic_algo': 'TNhaplotyper2'
        }
    
    # Set the job for production variant caller
    job_input['mappings_bam']= dx.dxlink(sample['sentieon-bwa_v3.2.0'])
    prod_job = dx.DXApp(name="sentieon-tnbam", alias="3.2.0").run(
        job_input,
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        folder="/sentieon-tnbam_v3.2.0",
        name=f"sentieon-tnbam_v3.2.0_{sample['Name']}"
        )
    job_ids['v3.2.0'] = prod_job.describe(fields={'id': True})['id']
    
    # Set the job for development variant caller
    job_input['mappings_bam']= dx.dxlink(sample['sentieon-bwa_v4.2.2'])
    dev_job = dx.DXApp(name="sentieon-tnbam", alias="5.0.1").run(
        job_input,
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        folder="/sentieon-tnbam_v5.0.1",
        name=f"sentieon-tnbam_v5.0.1_{sample['Name']}"
        )
    job_ids['v5.0.1'] = dev_job.describe(fields={'id': True})['id']
    
    # store job_ids in a list
    set_job_ids.append(job_ids)


json_jobs = json.dumps(set_job_ids)
print(json_jobs)