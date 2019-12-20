#!/bin/bash
# clone the directory
git clone https://github.com/nlapier2/Metalign.git
cd Metalign

# Download the data (this may take a few hours)
./setup_data.sh

# Create a conda virtual environment in order to install the tool locally
module load python/3.6.3-anaconda5.0.1  # load conda
conda create -y -n MetalignVE python=3.8  # create the environment
source activate MetalignVE  # activate the environment

# make sure you have access to GCC
module load gcc/7.3.1

# Install the Metalign requirements
./setup_libraries.sh
