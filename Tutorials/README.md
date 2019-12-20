# Tutorial
This directory contains an assortment of bash scripts that demonstrate (on a small amount of data) how to use this tools and a few different use cases for the pipeline.

## Local_install.sh
This shows how to perform a local installation of the tool.

## ICS_ACI_install.sh
This script shows how to perform a local installation on the Penn State ACI supercomputing cluster. It may be helpful for others running the code locally on a cluster (that does not allow/have docker).

## Run_full_pipeline.sh
A few different examples of how to run the full pipeline from the command line.

## Make_custom_database.sh
While we provide a pre-built database for you to use, you may be interested in utilizing a custom database. This script demonstrates the steps required to do this.

## Only_create_subset_database.sh
While we utilize Minimap2, you might want to use an alternate aligner. This script demonstrates how to select the relevant subset database using CMash, but not run the minimap2 step. The result will be a FASTA file that can be used as an alignment database for a different aligner.

## Map_but_dont_profile.sh
You may be interested in just the alignment and don't care about the taxonomic profile (or do not know or wish to create the CAMI formatted db_info.txt taxonomy file). This script shows how to run the method and output an alignment SAM file (without requiring the taxonomy information and without returning a taxonomic profile).

## Change_kmer_size.sh
The default analysis uses 60-mers: a reference organism will not be included in the subset database if it shares no 60-mer in common with the input sample data. You may wish to change this depending on the data you use (though our analyses show that k=60 is best for most data). This script demonstrates how to accomplish that.




