# A guide to searching for Relevant seqs

1. Pull all non-redundant assemblies from NCBI within the last 10 years
```
bash get-taxon-genomes.sh <NCBI_taxID>--include genome --mag exclude --exclude-multi-isolate --assembly-version latest --released-after <YYYY-MM-DD> --exclude-atypical
```