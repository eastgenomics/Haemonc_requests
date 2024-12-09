import dxpy
import json

cell_line_fastqc_list = [
    {"Name":"oncospan-cell-line-1st",
     "fastq_L001_R1":"file-G2X22Jj4qyfqzpybKfYJ15VQ",
     "fastq_L002_R1":"file-G2X25Q04qyfz1qf57gjxqx5f",
     "fastq_L001_R2":"file-G2X22KQ4qyfgpvyy9GQBkgF9",
     "fastq_L002_R2":"file-G2X25Qj4qyfZ77jf5kkYV88v"},
    {"Name":"oncospan-cell-line-2nd",
     "fastq_L001_R1":"file-G2X22PQ4qyfQZ9VvGQ5z2gyx",
     "fastq_L002_R1":"file-G2X25VQ4qyfV72xZJ63k9bJ4",
     "fastq_L001_R2":"file-G2X22Q04qyfqzpybKfYJ15Vz",
     "fastq_L002_R2":"file-G2X25XQ4qyfQZ9VvGQ5z2jF1"},
    {"Name":"HD734-1st",
     "fastq_L001_R1":"file-G2X22684qyfj1jjq60J56Z2V",
     "fastq_L002_R1":"file-G2X25984qyfZ399J3FB8Jgyp",
     "fastq_L001_R2":"file-G2X22704qyffpZXJ7vB5VVxX",
     "fastq_L002_R2":"file-G2X25B84qyfj9jbv39KFvYjy"},
    {"Name":"HD734-2nd",
     "fastq_L001_R1":"file-G2X228j4qyfkJ7zp28jVFYxV",
     "fastq_L002_R1":"file-G2X25F04qyff298GP8YPpbfg",
     "fastq_L001_R2":"file-G2X229Q4qyfgK0gJ8gb2v32K",
     "fastq_L002_R2":"file-G2X25G04qyfVvGYZ65YZyXQV"},
    {"Name":"HD829-1st",
     "fastq_L001_R1":"file-G2X222j4qyfvf2G21Z624ZYZ",
     "fastq_L002_R1":"file-G2X255Q4qyfQFGxJ21zy7X3j",
     "fastq_L001_R2":"file-G2X223Q4qyfz7VZB5F4yXbqG",
     "fastq_L002_R2":"file-G2X25684qyfb2FxpGfjg31x8"},
    {"Name":"HD829-2nd",
     "fastq_L001_R1":"file-G2X224Q4qyfVvGYZ65YZyX6k",
     "fastq_L002_R1":"file-G2X25784qyfvKQqY77Zj1pP6",
     "fastq_L001_R2":"file-G2X22584qyfpfZqq0PVBVkjZ",
     "fastq_L002_R2":"file-G2X25884qyfb7VV80gq5Xbby"}
    ]

sample_fastqc_list = [
    {"Name":"129404211-24107K0073-24NGSHO18-8128-M-96527893",
     "fastq_L001_R1":"file-Gjkp1B040738bfgj31Q6ZbJZ",
     "fastq_L002_R1":"file-Gjkp2y840739Qk7x7jX1Xjyj",
     "fastq_L001_R2":"file-Gjkp1Bj40734kp0qx349F0vP",
     "fastq_L002_R2":"file-Gjkp2zj40734j0ZP3ZJkbFz7"},
    {"Name":"129459591-24109K0055-24NGSHO18-8128-F-96527893",
     "fastq_L001_R1":"file-Gjkp18840738bfgj31Q6ZbJ0",
     "fastq_L002_R1":"file-Gjkp2pj40731vP0yzZ0Vq6yp",
     "fastq_L001_R2":"file-Gjkp18j40733v6f7z60y3PZX",
     "fastq_L002_R2":"file-Gjkp2q840733zj1x41yfjq8Y"},
    {"Name":"129450697-24109K0025-24NGSHO18-8128-M-96527893",
     "fastq_L001_R1":"file-Gjkp11j4073BFxbFy065B21X",
     "fastq_L002_R1":"file-Gjkp2P04073F222QQb5fb6Yb",
     "fastq_L001_R2":"file-Gjkp0yQ40734jgYjFp94Qx45",
     "fastq_L002_R2":"file-Gjkp2Pj407387fYz0973033Y"},
    {"Name":"129325254-24102K0083-24NGSHO17-8128-F-96527893",
     "fastq_L001_R1":"file-GjfY0Z04zZBGJG6G2972v2FZ",
     "fastq_L002_R1":"file-GjfY24j4zZB3719QZxf20gg4",
     "fastq_L001_R2":"file-GjfY0Z04zZBG4bZjXyKf75fB",
     "fastq_L002_R2":"file-GjfY24j4zZB7jzgYv35yy13j"}
]

all_fastqc_list = cell_line_fastqc_list + sample_fastqc_list


# Run all jobs and store all output jobs
set_job_ids = []
for sample in all_fastqc_list:
    job_ids = {}
    job_ids['Name'] = sample['Name']
    
    job_input = {
        "reads_fastqgzs": [dxpy.dxlink(sample['fastq_L001_R1']), dxpy.dxlink(sample['fastq_L002_R1'])],
        "reads2_fastqgzs": [dxpy.dxlink(sample['fastq_L001_R2']), dxpy.dxlink(sample['fastq_L002_R2'])],
        "genomebwaindex_targz": dxpy.dxlink("file-Gb76f204XGybZ3J6F731xkBp"),
        "genome_fastagz": dxpy.dxlink("file-Gb757784XGyY3FPvkPQ74K9z")
        }
    
    prod_job = dxpy.DXApp(name="sentieon-bwa", alias="3.2.0").run(
        job_input,
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        folder="/sentieon-bwa_v3.2.0",
        name=f"sentieon-bwa_v3.2.0_{sample['Name']}"
        )
    job_ids['v3.2.0'] = prod_job.describe(fields={'id': True})['id']
    
    dev_job = dxpy.DXApp(name="sentieon-bwa", alias="4.2.2").run(
        job_input,
        project="project-Gjz0P404fV6bq4J30qBfPBFf",
        folder="/sentieon-bwa_v4.2.2",
        name=f"sentieon-bwa_v4.2.2_{sample['Name']}"
        )
    job_ids['v4.2.2'] = dev_job.describe(fields={'id': True})['id']
    
    # store job_ids in a list
    set_job_ids.append(job_ids)
  

json_jobs = json.dumps(set_job_ids)
print(json_jobs)
