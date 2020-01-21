#!/bin/bash
# This script gives an example of how to install Metalign on a computing cluster via conda (i.e. when Docker is not available and you can't simply 'pip install' things). You may need to modify it depending on your cluster setup.

# After this install script, you need to be sure to activate the conda environment before running Metalign: source activate MetalignVE

# clone the directory
git clone https://github.com/nlapier2/Metalign.git
cd Metalign

# Create a conda virtual environment in order to install the CMash dependancy locally
git clone https://github.com/dkoslicki/CMash.git && cd CMash 
module load python/3.6.3-anaconda5.0.1  # load conda
conda env create -y -n MetalignVE -f tests/CMash-env-RedHat-Server-6.10.yml  # create the environment
source activate MetalignVE  # activate the environment

# you can optionally test if the CMash dependancy is installed correctly by running CMash/tests/./run_tests.sh

# make sure you have access to GCC
module load gcc/7.3.1

# Install the Metalign requirements
# ./setup_libraries.sh won't work in an environment where pip install can't/shouldn't be used or fails due to architecture, etc.
git clone https://github.com/lh3/minimap2.git && cd minimap2 && make && cd ..
git clone https://github.com/refresh-bio/KMC.git  && cd KMC && make && cd ..

# Download the data (this may take a few hours)
./setup_data.sh

