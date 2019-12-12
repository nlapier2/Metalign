# Data information

## spec_db_info.txt
This file describes the format of the required db_info.txt file. Note that a default db_info.txt file will be downloaded by running the `setup_data.sh` script.

## Default data after running the setup_data.sh script

An assortment of files are placed in this directory after running the `setup_data.sh` script. A brief description of each is as follows:

1. organism_files folder
This contains a database of 199,807 microbial organisms obtained from NCBI. The distribution of organisms is as follows:

| Kingdom    | Number of organisms in database |
|------------|---------------------------------|
|Bacteria    |                                 |
|Archaea     |                                 |
|Eukaryotes  |                                 |
|Viruses     |                                 |
|            |                                 |
|            |                                 |

2. db_info.txt
This file contains a CAMI formated taxonomy description of over 49M unique accessions obtained from every NCBI genome assembly as of November 2018.

3. cmash_db_n1000_k60.h5
This is an HDF5 formated file containing 1000 randomly selected (via MinHash) 60-mers from each of the organisms contained in the `organism_files` folder. Please see the [CMash repository](https://github.com/dkoslicki/CMash) for more information about how this was created (specifically, using the `CreateStreamingDNADatabase.py`).

4. cmash_db_n1000_k60.tst
This is a ternary search trie created from the randomly selected 60-mers from the training data. This was created using the [MakeStreamingDNADatabase.py](https://github.com/dkoslicki/CMash/blob/master/scripts/MakeStreamingDNADatabase.py) script in the CMash repository.

5. cmash_db_n1000_k60_dump.fa
This is a (FASTA formated) dump of all the selected 60-mers from the training data. This is used as input to KMC.

6. cmash_db_n1000_k60_dump.kmc_pre and cmash_db_n1000_k60_dump.kmc_suf
These are auxiliary files created by KMC when creating a k-mer search database from the cmash_db_n1000_k60_dump.fa file

7. cmash_filter_n1000_k60_30-60-10.bf and cmash_filter_n1000_k60_30-60-10.bf.desc
These are [hydra](https://github.com/crankycoder/hydra) bloom filters that are used as a quick pre-filter: we don't bother considering k-mers if they don't show up in the training data, and this filter identfies which those are.

