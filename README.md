# MiCoP: Microbial Community Profiling method for detecting viral and fungal organisms in metagenomic samples

MiCoP is a method for high-accuracy profiling of viral and fungal metagenomic communities. MiCoP uses a mapping-based approach that produces more accurate results relative to state of the art general purpose metagenomic profiling methods, which tend to focus mostly on bacteria. MiCoP is compatible with both single and paired end reads. This repository is set up to use BWA to map reads to the full NCBI Virus and Fungi RefSeq databases, and then filter and profile these results using the compute-abundances.py script. Usage is simple and details are provided below. For more details on MiCoP, see the paper link below.

### Synopsis

```
git clone https://github.com/smangul1/MiCoP.git
cd MiCoP/
./setup.sh
python run-bwa.py reads.fq [--virus OR --fungi] --output alignments.sam
python compute-abundances.py alignments.sam [--virus OR --fungi]
```

Where [--virus OR --fungi] means you use either the --virus flag or the --fungi flag (but not both or neither) depending on whether you want to profile viruses or fungi.

### Paper

The MiCoP paper is currently available on bioRxiv: https://www.biorxiv.org/content/early/2018/01/04/243188

If you use MiCoP, please cite the paper. For instance:

LaPierre, N., Mangul, S., Alser, M., Mandric, I., Wu, N. C., Koslicki, D., & Eskin, E. (2018). MiCoP: Microbial Community Profiling method for detecting viral and fungal organisms in metagenomic samples. *bioRxiv*, 243188.

### Download and Setup MiCoP

Simply run:
```
git clone https://github.com/smangul1/MiCoP.git
cd MiCoP/
./setup.sh
```

This will download this repository, change into the directory, and run the setup script. The setup script will download about 8.4 GB of data, which will be uncompressed into files totalling about 13 GB of disk space. These are the pre-indexed NCBI Virus and Fungi databases for BWA. It is not entirely necessary to use these specific databases, but it will greatly simplify running MiCoP and is strongly recommended.

This step should take about 5-30 minutes, depending on your network connection.

### Mapping stage

MiCoP is a mapping-based method that expects output from a mapping tool in SAM format. For the mapping tool, we recommend BWA. To eliminate confusion over which databases to use, how to extract necessary information from them, which alignment tool to use, and how to use it, we have included a copy of the BWA alignment tool (see License info) and a script to run it using the settings expected by MiCoP. This script is called run-bwa.py.

*Note*: Runtime for the BWA alignment step dominates the total computation time. For large reads files, this may take a while.

Basic usage is very simple:

```
python run-bwa.py reads.fq [--virus OR --fungi] --output alignments.sam
```

The reads are required, as are either --virus or --fungi, while the --output flag is optional. The default output file name is alignments.sam. If you are using paired end reads, use the same format as above except with two reads files ("python run-bwa.py reads1.fq reads2.fq" etc). If you have a single interleaved paired ends file, use the above format and also add the "--paired" flag.

If you would like to run BWA manually for more flexibility, make sure you use the -a flag, and see the following manpage for BWA: http://bio-bwa.sourceforge.net/bwa.shtml#13

Details on the basic BWA algorithm can be found in the BWA paper cited below. We specifically use the BWA mem algorithm, cited below BWA.

Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754-1760.

Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997.

### Abundance profiling Script

The script for computing abundance profiling is called compute-abundance.py, and is very simple to use:

```
python compute-abundances.py alignments.sam [--virus OR --fungi] [options]
```

This script should take a small fraction of the time it takes to run BWA mapping.

By default, this will output the results to abundances.txt. To change this, use the --output flag. All of the other options have to do with adjusting the filtering criteria for counting a clade as being present, and these options and their descriptions can be viewed with:

```
python compute-abundances.py -h
```

### Simulated data

The data used for our simulations can be found in the simulated\_data directory in a compressed format (.bz2). If you want to run BWA or the run-bwa.py wrapper script on this data, you will have to decompress it first. For instance, running "bzip2 -dk grinder-fungi-low-complexity-reads.fa.bz2" will generate the file grinder-fungi-low-complexity-reads.fa, which you can then run the scripts on. Just do "bzip2 -d" without the "k" if you only want the unpacked file and do not want to keep the original compressed file.

The mock community datasets are taken from the following papers:

Viral data: Conceicao-Neto, N., Zeller, M., Lefrere, H., De Bruyn, P., Beller, L., Deboutte, W., ... & Matthijnssens, J (2015). Modular approach to customise sample preparation procedures for viral metagenomics: a reproducible protocol for virome analysis. *Scientific reports*, 5, e16532. doi:10.1038/srep16532

Fungi data: Tonge, D. P., Pashley, C. H., & Gant, T. W. (2014). Amplicon–Based Metagenomic Analysis of Mixed Fungal Samples Using Proton Release Amplicon Sequencing. *PloS one*, 9(4), e93849. https://doi.org/10.1371/journal.pone.0093849

Finally, the HMP data was downloaded using the instructions in the MetaPhlAn tutorial: https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial

### License Info

For convenience, we include a version of the BWA kit in our repository. This is allowed under BWA's GPLv3 license, under the requirement that we also adopt GPLv3. We also use a code snippet in the setup script that is adapted from a Stack Overflow answer, permitted for inclusion under the Creative Commons license.
#
