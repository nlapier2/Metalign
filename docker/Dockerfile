# Define parent image
FROM ubuntu:18.04

# Default/home directory
WORKDIR /home

# Install needed packages
ENV DEBIAN_FRONTEND noninteractive
#RUN rm -rf /var/lib/apt/lists/*
RUN apt-get -qq update
RUN apt-get install -y curl build-essential git libz-dev vim zip
RUN apt-get install -y python-dev python-setuptools python-pip python-tk python3-dev
RUN python -m pip install numpy scipy virtualenv
RUN git clone https://github.com/dkoslicki/CMash.git && cd CMash && pip install -r requirements.txt && cd ..

# Copy current directory (including Metalign repo) into /home
ADD . /home

ENTRYPOINT ["python3"]
#
