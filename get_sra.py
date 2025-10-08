import csv
from Bio import Entrez
import sys
from argparse import ArgumentParser
from datetime import date
 # Example output: 2023-10-05

parser = ArgumentParser()
parser.add_argument("--email", dest="email", help="your.email@something.com as required by Entrez")
parser.add_Argument("--db", dest="organism", help="In the format <Genus_species>")

args = parser.parse_args() 
Entrez.email = args.email
taxa = args.organism
today = str(date.today()).replace("-","")


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
        #print(line)
        if any(tag in line for tag in ("SRR", "DRR", "ERR")):#if line.startswith("SRR"):
            srrs.append(line.split(",")[0])
    return list(set(srrs))


with open("passed_dist_filter.csv") as f:
    reader = csv.reader(f)
    with open(f"BigBacter.{taxa}.NCBI_filtered_{today}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample","taxa","genbank","sra"])
        for row in reader:
            filename, mash = row
            gcf_id = "_".join(filename.split("_")[:2])  # e.g. GCF_002120665.1
            biosample = get_biosample_from_gcf(gcf_id)
            if biosample:
                srrs = get_srr_from_biosample(biosample)
                if srrs:
                    for srr in srrs:
                        print(gcf_id, biosample, srr, mash)
                        writer.writerow([biosample, taxa, gcf_id, srr])   # wrap in list so each line is one cell
                else: pass#print(gcf_id,biosample, srrs)
            else:
                print(gcf_id, "No BioSample found")

