#/usr/bin/env xonsh
"""
Based on code examples here:
    https://github.com/schultzm/entrez_direct_tut
"""

from pathlib import Path
from collections import defaultdict
import pandas as pd
import click
import csv
import os

@click.command()
@click.option('-p', '--bioproject', 'bioproject', required=True, type=str, help = "BioProject accession")
@click.option('-f', '--outdirfiles', 'outdirfiles',default="files", show_default=True, type=str, help = "Name of directory to save files")
@click.option('-s', '--outdirsras', 'outdirsras',default="srafiles", show_default=True , type=str, help = "Name of directory to save .sra files")
@click.option('-fq', '--outdirfastqs', 'outdirfastqs',default="fastqs", show_default=True , type=str, help = "Name of directory to save .fastq files")
@click.option('-l', '--outdirlogs', 'outdirlogs',default="logs", show_default=True , type=str, help = "Name of directory to save log files")

def getsra(bioproject, outdirfiles, outdirsras, outdirfastqs, outdirlogs):
    workdir= os.getcwd()
    sample_file = outdirfiles + "/" + bioproject + "_biosamples.txt"
    srs_file = outdirfiles + "/" + bioproject + "_SRStoSRR.txt"
    log = workdir + "/" + outdirlogs + "/prefetched.log"
    os.makedirs(os.path.dirname(log), exist_ok=True)
    touch @(log)
    os.makedirs(os.path.dirname(outdirfiles + '/'), exist_ok=True)
    os.makedirs(os.path.dirname(outdirsras + '/'), exist_ok=True)
    os.makedirs(os.path.dirname(outdirfastqs + '/'), exist_ok=True)

    # Get Accession numbers, Names, and SRS numbers for each sample associated with the bioproject and write this info to a file
    print("Getting Sample Info")
    if not Path(sample_file).exists():
        esearch -db bioproject -query @(bioproject) | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -block Ids -element Id -group sra > @(sample_file)
    else:
        print("Skipping step, " + sample_file + " already exists")

    # Read the Sample file as a data frame and filter out samples without SRS accession
    df = pd.read_table(sample_file, header=None, names=["SAM", "Name", "SRS"])
    df.dropna(subset=['SRS'], inplace=True)

    # Get sequencing run IDs associated with each sample (Don't reconstruct SRS to SRR table if it already exists)
    dsrs = defaultdict(list)

    print("Getting Run Info")
    if not Path(srs_file).exists():
        ofile = open(srs_file,"w")
        for srs in df.SRS:
            print(f"Processing {srs}")
            srr = $(esearch -db sra -query @(srs) | efetch -format runinfo -mode xml | xtract -pattern Run -element Run)
            for item in srr.split("\n"):
                if not len(item): 
                    continue
                dsrs[srs].append(item)
                ofile.write(f"{srs}\t{item}\n")
            ofile.flush()
        ofile.close()
    else:
        dfsrs = pd.read_table(srs_file, header=None, names = ['SRS', 'SRR'])
        for (srs,srr) in zip(dfsrs.SRS, dfsrs.SRR):
            dsrs[srs].append(srr)

    # Prefetch the SRA files
    print("Start Prefetching")
    for sample, run in dsrs.items():
        print(f"Prefetching {run}")
        mkdir -p @(outdirsras + "/" + sample)
        cd @(outdirsras + "/" + sample)
        prefetch @(run) &>> @(log)
        cd ../..

	# Download again if run was not successfully downloaded
    print("Start Prefetching of failed samples")
    max_iterations = 5
    iterations = 0
    failed = ['initiate loop']
    while failed and iterations < max_iterations:
        success = $(grep "downloaded successfully" @(log) | cut -d " " -f4).split('\n')
        success = [s.replace("'", "") for s in success]
        success = [item for item in success if item != ""]
        all_values = [value for sublist in dsrs.values() for value in sublist]
        failed = [value for value in all_values if value not in success]
        failed_dsrs = {key: value for key, value_list in dsrs.items() for value in value_list if value in failed}

        for sample, run in failed_dsrs.items():
                print(f"Prefetching again : {run}")
                mkdir -p @(outdirsras + "/" + sample)
                cd @(outdirsras + "/" + sample)
                prefetch @(run) &>> @(log)
                cd ../..
        iterations += 1

    # Fasterq-dump
    for each in Path(outdirsras).resolve().rglob("*.sra"):
        ck = Path(outdirfastqs) / f"{each.parent.name}_1.fastq"
        if ck.exists(): 
            continue
        print(f"Processing {each}")
        fasterq-dump @(str(each)) -O @(outdirfastqs) -p -S

    # Getting read pair tables
    print("Getting read pair tables")
    paired = []
    unpaired = []

    for srs in sorted(Path(outdirsras).glob("SRS*")):
        for srr in sorted(srs.glob("SRR*")):
            r1 = Path(outdirfastqs) / f"{srr.name}_1.fastq"
            r2 = Path(outdirfastqs) / f"{srr.name}_2.fastq"
            if r1.exists() and r2.exists(): 
                size1 = r1.stat().st_size
                size2 = r2.stat().st_size
                paired.append((srs.name, srr.name, r1.name, r2.name, size1+size2))
            else:
                ufiles = ",".join([str(i.name) for i in Path(outdirfastqs).glob(f"{srr.name}*.fastq")])
                unpaired.append((srs.name, srr.name, ufiles))

    with open(outdirfiles + '/read_pair_table.csv', 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerow(("sample", "run", "file1", "file2", "size"))
        w.writerows(paired)


    with open(outdirfiles + '/unpaired_fastqs.csv', 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerow(("sample", "run", "files"))
        w.writerows(unpaired)

    SRS = [i[0] for i in paired]
    setSRS = set(SRS)

    with open(outdirfiles + "/samples.txt", "w") as outfile:
        outfile.write("\n".join(setSRS))

if __name__ == "__main__":
    getsra() # Runs the function
