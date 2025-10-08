# A guide to searching for Relevant seqs to process with BigBacter
This process is designed with BigBacter usage in mind, but scripts can be tweaked to work for other analyses. 

1. Pull all non-redundant assemblies from NCBI within the last 10 years
```
bash get-taxon-genomes.sh <NCBI_taxID>--include genome --mag exclude --exclude-multi-isolate --assembly-version latest --released-after <YYYY-MM-DD> --exclude-atypical
```


2. Get the relevant cluster reference genomes. We are going to compare the downloaded NCBI assemblies from NCBI
```
python get_bb_references.py --bucket <s3_bucket> --db <Genus_species>
```

This will pull down BigBacter cluster reference genomes will be saved as `bb_cluster_references/<Genus_species>_<cluster#>_ref.fa.gz`


3. Make a reference mash database from the cluster reference genomes
```
bash make_ref_db.sh
```
This will output a ref mash db as `ref_mash_db.msh`


4. Determine mash distances between each downloaded NCBI assembly and each cluster reference genome. First, you must determine the mash threshold to use for your organism. We recommend setting the scope of your filtering to still keep sequences that might not be close enough related to cluster the investigation sequences (=<100 SNPs) for two reasons; to provide context of the mash distances within your investigative clusters vs external NCBI sequences and to give buffer for mash distance due to recombination which BigBacter can mask when needed.

A mash distance of 0.001 - 0.005 should work for most gram positive bacteria. Verify that `mash_distance*genome_size` is greater than 100 SNPs for your organism. 

***NOTE*** This is a starting estimate mash distance, increase the allowed distance if no sequences pass your distance filter. Additional consideration should be taken for organisms that are known to be highly recombinate such as Legionella.
For example with Staphylococcus_aureus

```
# Staphylococcus_aureus genome size = 2,800,000
0.001*2,800,000 = 2,800 # Greater than 100
```

This will output two files of interest: 

-`ref_nearest_nonref.csv` <br>
    Contains the mash distance to for each cluster reference genome to the NCBI assembly that is most closely related. This helps provide context of how significant the within cluster SNP distances are to the distances to sequences outside of the cluster.

-`passed_dist_filter.csv` <br>
    CSV that contains the GenBank and/or RefSeq IDs that are equal to or less than the chosen mash distance to any of the cluster reference genomes. The mash distance listed for each ID is the smallest mash distance determined for the NCBI genome and any of the cluster reference genomes.

After completing this step, view the contexts of the files generated. If no NCBI sequences are output to `passed_dist_filter`, you may need to increase the allowed mash distance or reducing the filtering in step 1 (and pull another dataset).


5. BigBacter requires both GenBank/RefSeq ID and SRR as input when using the `--ncbi` arg. Previous steps have been working soley with the GenBank/RefSeq ID, next is pulling the associated SRR _if one exists_. Many assemblies on NCBI *do not* have associated SRR, thus missing SRRs is not a bug but merely reflective to what actually exists. It is possible that no SRRs get returned in which case the data for possibly closely related seqs is not publically available on NCBI. To determine the (possible) SRR IDs for the NCBI assemblies, run:

```
python get_sra.py --email <your_email@something.com> --db <BB_db>
```

The output of this process will be `BigBacter.<db>.NCBI_filtered_<today_YYYYMMDD>.csv`