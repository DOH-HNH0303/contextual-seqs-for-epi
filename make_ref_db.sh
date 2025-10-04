#!/bin/bash

rm -rf bb_ref_list.txt
rm ref_mash_db*
ls bb_cluster_references/ | sed 's|^|bb_cluster_references/|' > bb_ref_list.txt
mash sketch -o ref_mash_db -l bb_ref_list.txt
