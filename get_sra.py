import csv
from Bio import Entrez
import sys
import argparse
from argparse import ArgumentParser
from datetime import date
import json
 # Example output: 2023-10-05

parser = ArgumentParser()
parser.add_argument("--email", dest="email", help="your.email@something.com as required by Entrez")
parser.add_argument("--db", dest="organism", help="In the format <Genus_species>")
parser.add_argument('--illumina', dest='illumina', action=argparse.BooleanOptionalAction)
parser.add_argument('--ont', action=argparse.BooleanOptionalAction)
parser.add_argument('--pacbio', action=argparse.BooleanOptionalAction)
# TODO  -add args for --ont,--illumina,--pacbio
args = parser.parse_args() 
Entrez.email = args.email
taxa = args.organism
illumina = args.illumina
ont = args.ont
pacbio = args.pacbio

tech_list = []
if illumina:
    tech_list.append("illumina")
if pacbio:
    tech_list.append("pacbio")
if ont:
    tech_list.append("ont")


today = str(date.today()).replace("-","")
assembly_report_file="assemblies/assembly_data_report.jsonl"
with open(assembly_report_file, 'r') as f:
            for line in f:
                record = json.loads(line)


def check_for_any_substring(target_string, substring_list):
    """
    Checks if any substring from a list is present in a target string.

    Args:
        target_string (str): The string to search within.
        substring_list (list): A list of substrings to check for.

    Returns:
        bool: True if any substring from the list is found in the target_string,
              False otherwise.
    """
    return any(sub in target_string for sub in substring_list)


def get_seq_technology(report_file=assembly_report_file,tech_list=tech_list):
    """
    Filters a list of NCBI accession IDs based on the sequencing technology used
    for the assembly, using the assembly_data_report.jsonl file.

    Args:
        input_ids_file (str): Path to a file containing one accession ID per line.
        assembly_report_file (str): Path to the assembly_data_report.jsonl file.
        output_pacbio_file (str): Path to output file for PacBio IDs.
        output_ont_file (str): Path to output file for ONT (Oxford Nanopore) IDs.
        output_illumina_file (str): Path to output file for Illumina IDs.
    """

    # Dictionary to store the sequencing technology for each assembly
    assembly_tech_map = {}
    # tech_list = []
    # if illumina:
    #     tech_list.append("illumina")
    # if pacbio:
    #     tech_list.append("pacbio")
    # if ont:
    #     tech_list.append("ont")
    
    # Parse the assembly report to build the technology map
    try:
        with open(assembly_report_file, 'r') as f:
            for line in f:
                record = json.loads(line)
                # if "GCF_900097465.1" in str(record):
                #     print(record)
                if 'sequencingTech' in record['assemblyInfo'].keys() or 'sequencingTech' in record.keys():
                    
                    if 'sequencingTech' in record.keys():
                        tech = record['sequencingTech'].split("; ")
                    else:
                        tech = record['assemblyInfo']['sequencingTech'].split("; ")
                    tech = [w.lower() for w in tech]
                    tech_clean = []
                    for s in tech:
                        # collect only allowed substrings that appear in the string
                        keep = [a for a in tech_list if a in s.lower()]
                        # join them back into a string (or empty if none matched)
                        tech_clean.append(" ".join(keep))
                        
                    if set(tech_list) & set(tech_clean):
                        genbank = record['accession']
                        if isinstance(assembly_tech_map, dict):
                            if genbank not in assembly_tech_map:
                                assembly_tech_map[record['accession']] = [tech_clean, record['assemblyInfo']['biosample']['accession']]
                        else:
                            assembly_tech_map[record['accession']] = [tech_clean, record['assemblyInfo']['biosample']['accession']]
                        if record['pairedAccession']:
                            if assembly_tech_map:
                                if record['pairedAccession'] not in assembly_tech_map:
                                    assembly_tech_map[record['pairedAccession']] = [tech_clean, record['assemblyInfo']['biosample']['accession']]
                            else:
                                assembly_tech_map[record['pairedAccession']] = [tech_clean, record['assemblyInfo']['biosample']['accession']]

    except FileNotFoundError:
        print(f"Error: Assembly report file not found at '{assembly_report_file}'")
        return
    except json.JSONDecodeError as e:
        print(f"Error: Could not parse JSON in '{assembly_report_file}'. {e}")
        return
    return assembly_tech_map



def get_biosample_from_gcf(gcf_id):
    """Return the BioSample accession (SAMN...) for a given assembly accession."""
    handle = Entrez.esearch(db="assembly", term=gcf_id)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return None
    uid = record["IdList"][0]

    # Fetch the assembly summary
    summary = Entrez.esummary(db="assembly", id=uid)
    docsum = Entrez.read(summary)
    summary.close()
    #print(docsum['DocumentSummarySet']['DocumentSummary'][0].get("BioSampleAccn"))
    biosample = docsum['DocumentSummarySet']['DocumentSummary'][0].get("BioSampleAccn")
    return biosample

def get_srr_from_biosample(biosample_accn):
    """Return a list of SRR IDs linked to a BioSample accession."""
    srrs = []
    # Search SRA for this BioSample
    handle = Entrez.esearch(db="sra", term=biosample_accn)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return srrs

    # Fetch runinfo for all linked SRA records
    ids = ",".join(record["IdList"])
    runinfo = Entrez.efetch(db="sra", id=ids, rettype="runinfo", retmode="text")
    lines = runinfo.read().splitlines()
    runinfo.close()

    #print(lines)
    for line in lines:
        line = line.decode("utf-8")
        if check_for_any_substring(str(line).lower(), tech_list):
            print()
            #print(line)
            if any(tag in line for tag in ("SRR", "DRR", "ERR")):#if line.startswith("SRR"):
                srrs.append(line.split(",")[0])
    return list(set(srrs))


with open("passed_dist_filter.csv") as f:
    reader = csv.reader(f)
    data_assembly_json = get_seq_technology()
    #exit(1)
    with open(f"BigBacter.{taxa}.NCBI_filtered_{today}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample","taxa","genbank","sra"])
        for row in reader:
            filename, mash = row
            if filename == "assembly_file":
                pass
            else:
                gcf_id = "_".join(filename.split("_")[:2])  # e.g. GCF_002120665.1
                if gcf_id in data_assembly_json:
                    print(data_assembly_json[gcf_id])
                    biosample = data_assembly_json[gcf_id][1]
                    print(biosample)
                    if biosample:
                        srrs = get_srr_from_biosample(biosample)
                        if srrs:
                            for srr in srrs:
                                print(gcf_id, biosample, srr, mash)
                                writer.writerow([biosample, taxa, gcf_id, srr])   # wrap in list so each line is one cell
                        else: pass#print(gcf_id,biosample, srrs)
                    else:
                        print(gcf_id, "No BioSample found")

