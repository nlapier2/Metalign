## Overview

This directory contains experimental results and a script to reproduce the figures used in the paper (paper_plots.ipynb).

The figures/ directory contains the figures used in the paper.

The raw_results/ directory contains the raw output results from the programs run.

The scripts/ directory contains the scripts submitted to the job cluster and other lines run to generate the raw results. It also has scripts to process these results into a format that is easier to make the figures from. Those reformatted files are in the processed_results/ directory. 


## Data download instructions

We exclusively used publicly-available datasets in this paper. For the CAMI 1 results, we used the results and datasets from the GigaDB dataset here: http://dx.doi.org/10.5524/100344. We used the challenge datasets, not the toy datasets. For the medium complexity datasets, we used the insert size 270 reads. 

The Peabody mock community data was retrieved using the link in their paper (see: https://doi.org/10.1186/s12859-015-0788-5). 

The prokaryote-isolated Tara Oceans reads used in this study are available on EBI here: https://www.ebi.ac.uk/ena/data/view/PRJEB1787. The run accessions we used, which were not chosen for any particular reason, were ERR598948, ERR598949, ERR598950, ERR598952, ERR598957. Specifically, we went to the first page in which the Instrument Model for all samples was Illumina HiSeq 2000. We then downloaded all ten datasets on that page, and decided to use only the deep chlorophyll maximum samples, so that the samples would be from a consistent ocean layer. Extra details on each file (sampling location, depth, isolation method, etc) can be found by clicking on the sample accession corresponding to each of those run accessions. Information on the locations of some of the Tara Oceans stations can be found here: http://ocean-microbiome.embl.de/companion.html.
