# Define parent image
FROM ubuntu:18.04

# Default/home directory
WORKDIR /home

# Install needed packages
ENV DEBIAN_FRONTEND noninteractive
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get -qq update
RUN apt-get install -y curl build-essential git libz-dev vim zip
RUN apt-get install -y python-dev python-setuptools python-pip python-tk python3-dev
RUN python -m pip install numpy scipy
RUN git clone https://github.com/nlapier2/Metalign.git && cd Metalign && ./setup_libraries.sh
# A quick way to update the dockerfile
RUN cd Metalign/ && git pull origin master && echo "Metalign docker version 0.4.0"

ENTRYPOINT ["python3"]
#
