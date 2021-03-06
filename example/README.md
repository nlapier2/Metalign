## Running some example data

We use a small example dataset that users can download and use to test that they have installed Metalign properly. This data is the mock community data used in the paper. The data was generated by Peabody et al. [1] and is in the public domain. Please see example/data_license_info for more information.

[1] Peabody, M.A., Van Rossum, T., Lo, R. et al. Evaluation of shotgun metagenomics sequence classification methods using in silico and in vitro simulated communities. BMC Bioinformatics 16, 362 (2015) doi:10.1186/s12859-015-0788-5

## Instructions

The instructions below assume your current directory is the example/ directory.

Download the data:

` wget https://ucla.box.com/shared/static/ybz1xgke32kh56p4lqsg41t4g1bq8xy0.gz && mv ybz1xgke32kh56p4lqsg41t4g1bq8xy0.gz reads.fna.gz `

Run metalign via docker:

` docker run --rm -it -v /path/to/Metalign/:/home/workdir/ nlapier2/metalign Metalign/metalign.py /home/workdir/example/reads.fna.gz /home/workdir/data/ --output /home/workdir/example/abundances.tsv `

OR

Run metalign via local repository:

` python3 ../metalign.py reads.fna.gz ../data/ `

This should take about 10-15 minutes. Check to make sure that the abundances.tsv file looks reasonable. If you run

` head abundances.tsv `

you should get results that look roughly like:

```
@SampleID:/home/workdir/example/reads.fna.gz
@Version:Metalign-v0.2
@Ranks: superkingdom|phylum|class|order|family|genus|species|strain

@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE	_CAMI_genomeID	_CAMI_OTU
2	superkingdom	2	Bacteria	99.65183
1224	phylum	2|1224	Bacteria|Proteobacteria	63.12215
201174	phylum	2|201174	Bacteria|Actinobacteria	24.41209
1239	phylum	2|1239	Bacteria|Firmicutes	12.11759
1236	class	2|1224|1236	Bacteria|Proteobacteria|Gammaproteobacteria	47.25683
```
