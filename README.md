# Metalign: Efficient alignment-based metagenomic profiling via containment min hash

Metalign is a method for estimating the taxonomic composition and relative abundances of organisms in a metagenomic sample based on whole-genome shotgun sequencing reads.

### Quick Usage

```
git clone https://github.com/nlapier2/Metalign.git
cd Metalign/
./setup.sh
python3 metalign.py my_reads.fq --output profile.tsv
```

Note that this will run some python package installations via "pip install --user". If you want to root install, run "./setup.sh root" instead. Alternatively, you can install the bioconda package or manually set up the package. These options are described below.

### Paper

The Metalign paper is currently available on bioRxiv: (Not yet)

If you use Metalign, please cite the paper. For instance:

LaPierre, N., Koslicki, D., Alser, M., Eskin, E., & Mangul, S. (2019). Metalign: Efficient alignment-based metagenomic profiling via containment min hash. *bioRxiv*, (Not yet).

### Download and Setup

Simply run:
```
git clone https://github.com/nlapier2/Metalign.git
cd Metalign/
./setup.sh
```

This will download this repository, change into the directory, and run the setup script.

This step should take a few minutes up to an hour or so, depending on your network connection and how many of the CMash dependenies you already have installed.

(Describe alternate steps here)

(More sections covering each step and detailed usage)
#
