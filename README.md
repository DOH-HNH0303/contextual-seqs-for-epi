# 🧬 Guide to Identifying Relevant Sequences for BigBacter Analysis

This workflow is designed for use with **BigBacter**, but the scripts can be adapted for other comparative genomics pipelines.

---

## 📥 Step 1: Download Non-Redundant Assemblies from NCBI

Use the following script to pull all non-redundant genome assemblies released in the last 10 years:

```
bash get-taxon-genomes.sh <NCBI_taxID> \
  --include genome \
  --mag exclude \
  --exclude-multi-isolate \
  --assembly-version latest \
  --released-after <YYYY-MM-DD> \
  --exclude-atypical
```

---

## 📦 Step 2: Retrieve Cluster Reference Genomes

Download BigBacter cluster reference genomes for your target genus/species:

```
python get_bb_references.py --bucket <s3_bucket> --db <Genus_species>
```

Output files will be saved to:

`bb_cluster_references/<Genus_species>_<cluster#>_ref.fa.gz`

---

## 🗃️ Step 3: Build Mash Reference Database

Create a Mash database from the reference genomes:

```
bash make_ref_db.sh
```

Output:

`ref_mash_db.msh`

---

## 📏 Step 4: Calculate Mash Distances

Compare downloaded NCBI assemblies to cluster reference genomes using Mash. First, determine an appropriate Mash distance threshold for your organism.

### 🔬 Recommended Thresholds

- For most Gram-positive bacteria: 0.001 - 0.005
- Ensure: mash_distance × genome_size > 100 SNPs

Example for *Staphylococcus aureus* (genome size ≈ 2.8 Mb):

0.001 × 2,800,000 = 2,800 SNPs

### ⚠️ Notes

- If no sequences pass the filter, increase the Mash distance or relax filters in Step 1.
- Highly recombinant organisms (e.g., *Legionella*) may require higher thresholds.

Run the `get_mash_dist.sh` script; use `tmux` if necessary:

```
bash get_mash_dist.sh <mash_distance threshold;default is 0.001>
```


### 📂 Output Files

- `ref_nearest_nonref.csv`  
  Mash distance between each reference genome and its closest NCBI match.

- `passed_dist_filter.csv`  
  GenBank/RefSeq IDs within the chosen Mash distance threshold.

---

## 🔍 Step 5: Retrieve SRR IDs for Filtered Assemblies

BigBacter requires both GenBank/RefSeq ID and SRR ID when using the `--ncbi` argument. Many assemblies lack associated SRRs—this is expected and not a bug.

To retrieve SRR IDs (if available):

```
python get_sra.py --email <your_email@domain.com> --db <BB_db>
```

Output:

`BigBacter.<db>.NCBI_filtered_<YYYYMMDD>.csv`

---

## 🙌 Acknowledgments

Special thanks to the contributors who made BigBacter possible:

|Name|Association|Contribution|
|-|-|-|
|[Shawn Hawkens](https://github.com/@DOH-SEH2303) |WA DOH, Molecular Epi.|Co-contributor for which this project would not exist without|
|[Jared Johnson](https://github.com/@DOH-JDJ0303) |WA DOH, Bioinformatics|Guidance with BigBacter and valueable sounding board for ideas|
|[Dahlia Walters](https://github.com/@DOH-DEW3303)|WA DOH, Molecular Epi.|Provides guidance to local scope for which this project is targeted|

If your name is missing or information is incorrect, please contact:  
waphl-bioinformatics@doh.wa.gov
